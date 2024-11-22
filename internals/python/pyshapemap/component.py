# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import subprocess as sp
import shutil
import signal
import os

# FIXME: this is just to avoid name collisions between connect() from connect.py
#        and connect() from nodes.py. Gotta be a better way
from pyshapemap.nodes import Node, InputNode, OutputNode, ParameterNode
from pyshapemap.nodes import FileNode, FolderNode, SharedInputNode, Edge
from pyshapemap.nodes import StdoutNode, StderrNode, StdinNode, ComponentNode
from pyshapemap.nodes import PipeNode

from pyshapemap.util import rand_id

# FIXME: rename most funcs to be "private" (i.e. __whatever__)
# FIXME: add __delattr__ and __delitem__ to Component class

class Component():
    def __init__(self, **kwargs):
        #####
        # These properties are initialized in overridden __setattr__()
        # to avoid problems with subclass constructors
        # self.internal_components = []
        # self.input_nodes = []
        # self.output_nodes = []
        #####
        self.name = None
        self.assoc_rna = None
        self.assoc_sample = None
        self.parent_component = None
        self.proc = None
        self.hung = False
        self.id = rand_id()
        self.progmon = None
        self.gv_props = {"style": '',
                         "color": 'black',
                         "shape": 'box',
                         "fontsize": '24'}
        if "internal_components" in kwargs:
            self.add(kwargs["internal_components"])
            del kwargs["internal_components"]
        if "gv_props" in kwargs:
            self.gv_props.update(kwargs["gv_props"])
            del kwargs["gv_props"]
        # WIP: eliminating silent parameter passthrough
        #      (resulted in bugs in the past if a parameter name
        #       changed and was not updated in some constructor)
        for arg in kwargs:
            if arg in self.__dict__:
                self.__dict__[arg] = kwargs[arg]
            else:
                raise RuntimeError("Unrecognized parameter \"{}\"".format(arg))
        
    def __getattr__(self, name):
        """
        Allow simple component.node_name access to nodes or components
        inside a component.

        """
        # Note: python should only be calling this function after
        #       checking if name is present in self.__dict__
        for node in self.input_nodes:
            if name == node.name:
                return node
        for node in self.output_nodes:
            if name == node.name:
                return node
        for comp in self.internal_components:
            if name == comp.get_name():
                return comp
        raise AttributeError

    def __getitem__(self,
                    name):
        """
        Also allow container-like access to internal components and nodes.
        For example: comp.MutationParser["parsed"]
        """
        try:
            assert isinstance(name, str)
        except AssertionError:
            pass
        return self.__getattr__(name)


    def __setattr__(self, name, value):
        found_name = False

        # prevent infinite recursions by making sure certain members are
        # initialized before attempting to access
        if "internal_components" not in self.__dict__:
            self.__dict__["internal_components"] = []
        if "input_nodes" not in self.__dict__:
            self.__dict__["input_nodes"] = []
        if "output_nodes" not in self.__dict__:
            self.__dict__["output_nodes"] = []

        for i in range(len(self.internal_components)):
            if name == self.internal_components[i].get_name():
                self.internal_components[i] = value
                value.parent_component = self # This is important for reorganizing existing pipelines
                found_name = True
        for i in range(len(self.input_nodes)):
            if name == self.input_nodes[i].name:
                self.input_nodes[i] = value
                # only set parent_component of the node if this is a
                # lowest-level component
                if len(list(self.internal_components)) == 0:
                    value.parent_component = self
                found_name = True
        for i in range(len(self.output_nodes)):
            if name == self.output_nodes[i].name:
                self.output_nodes[i] = value
                # only set parent_component of the node if this is a
                # lowest-level component
                if len(list(self.internal_components)) == 0:
                    value.parent_component = self
                found_name = True
        if not found_name:
            super().__setattr__(name, value)

    def __setitem__(self, name, value):
        """
        Also allow container-like access
        """
        assert isinstance(name, str)
        self.__setattr__(name, value)

    def get_name(self):
        name = self.__class__.__name__
        if self.name is not None:
            name = self.name
        return name

    def get_parent_names(self):
        # return list of all parent component names
        names = []
        p = self.parent_component
        while p is not None:
            names.append(p.get_name())
            p = p.parent_component
        return names


    def add(self,
            o,
            alias=None):
        """

        Args:
            o: InputNode/OutputNode/ParameterNode (to be appended to self.input_nodes
                     or self.output_nodes), Component (to be appended
                     to self.internal_components), or list

        """
        if isinstance(o, list):
            for x in o:
                self.add(x)
            return
        elif isinstance(o, Component):
            self.add_internal_component(o)
        elif isinstance(o, (InputNode, OutputNode, ParameterNode)):
            self.add_node(o, alias=alias)
        else:
            raise RuntimeError("Error: unrecognized object type \"{}\" for Component add()".format(type(o).__name__))

    def insert(self,
               to_insert,
               downstream_object):
        if isinstance(to_insert, list):
            for x in to_insert:
                self.insert(x, downstream_object)
            return
        assert isinstance(to_insert, Component)
        assert isinstance(downstream_object, Component)
        assert downstream_object in self.internal_components
        if to_insert in self.internal_components:
            raise RuntimeError("Error: component \"{}\" is already present in component \"{}\"".format(to_insert.get_name(), self.get_name()))
        if to_insert.get_name() in [c.get_name() for c in self.internal_components]:
            raise RuntimeError("Error: component \"{}\" has the same name as another component already present in component \"{}\". Give the component a unique name.".format(to_insert.get_name(),self.get_name()))
        o = to_insert
        d = downstream_object
        self.internal_components.insert(self.internal_components.index(d), o)
        o.parent_component = self

    def remove(self,
               *args):
        """
        Remove this Component from its parent (if no args) or
        remove provided object from internal_components, input_nodes,
        or output_nodes depending on object type.

        """
        assert len(args)==0 or len(args)==1
        if len(args)==0:
            self.remove_from_parent()
        elif len(args)==1:
            self.remove_internal_component(args[0])
        else:
            raise RuntimeError("Too many args provided to Component remove()")

    def replace(self,
                c1, c2):
        """
        Replace Component c1 with c2 in internal_components

        """
        assert isinstance(c1, Component)
        assert isinstance(c2, Component)
        try:
            index = self.internal_components.index(c1)
        except ValueError:
            raise RuntimeError("Component {} not present in {}".format(c1, self))
        self.internal_components[index] = c2

    def add_internal_component(self, component):
        assert isinstance(component, Component)
        if component in self.internal_components:
            raise RuntimeError("Error: component \"{}\" is already present in component \"{}\"".format(component.get_name(), self.get_name()))
        if component.get_name() in [c.get_name() for c in self.internal_components]:
            raise RuntimeError("Error: component \"{}\" has the same name as another component already present in component \"{}\". Components must have unique names.".format(component.get_name(),self.get_name()))
        self.internal_components.append(component)
        component.parent_component = self

    def remove_from_parent(self):
        """
        Remove this component from the internal_components list of its
        parent_component.
        """
        assert self.parent_component is not None
        assert self in self.parent_component.internal_components
        self.parent_component.internal_components.pop(self.parent_component.internal_components.index(self))
        self.parent_component = None

    def remove_internal_component(self,
               o):
        assert isinstance(o, Component)
        self.internal_components.pop(self.internal_components.index(o))
        o.parent_component = None

    def remove_node(self,
                         node):
        assert isinstance(node, (InputNode, OutputNode))
        if isinstance(node, InputNode):
            self.input_nodes.pop(self.input_nodes.index(InputNode))
        elif isinstance(node, OutputNode):
            self.output_nodes.pop(self.output_nodes.index(OutputNode))
        node.parent_component = None

    def add_node(self,
                 node,
                 alias=None):
        assert isinstance(node, (InputNode, OutputNode, ParameterNode))
        if isinstance(node, (InputNode, ParameterNode)):
            self.add_input_node(node, alias=alias)
        elif isinstance(node, OutputNode):
            self.add_output_node(node, alias=alias)


    def add_input_node(self, node, alias=None):
        assert isinstance(node, (InputNode, ParameterNode))
        self.input_nodes.append(node)
        if alias is not None:
            assert isinstance(alias, str)
            # allow access to this node using a different name
            # (useful if, for example, a component has a single
            # main output that may depend on different internal
            # components depending on the number of inputs)
            self.__dict__[alias] = node
        # only set parent_component of the node if this is a
        # lowest-level component and the node does not already
        # have a parent_component
        if ( len(self.internal_components) == 0 and
             node.parent_component is None ):
            node.parent_component = self

    def add_output_node(self, node, alias=None):
        assert isinstance(node, OutputNode)
        self.output_nodes.append(node)
        if alias is not None:
            assert isinstance(alias, str)
            # allow access to this node using a different name
            self.__dict__[alias] = node
        # only set parent_component of the node if this is a
        # lowest-level component and the node does not already
        # have a parent_component
        if ( len(self.internal_components) == 0 and
             node.parent_component is None ):
            node.parent_component = self

    def get_component_nodes(self):
        return self.input_nodes + self.output_nodes

    def format_command(self,
                       command):
        """
        Replace filename placeholders with quoted filenames

        """
        if command is None:
            msg = "Error: component \"{}\" cmd() does not return"
            msg += " a command (either string or list of strings)."
            msg = msg.format(self.get_name())
            raise RuntimeError(msg)


        # if command is given as a list of str args, combine into
        # a single string
        if isinstance(command, list):
            assert all((isinstance(x,str) for x in command))
            command = ' '.join(command)
        else:
            assert isinstance(command, str)

        values = {}
        for node in self.get_component_nodes():
            if node.get_name() in ["stdout", "stderr"]:
                continue
            # FIXME: all these isinstances are hacky, should use a unified interface to access file/folder path, move logic to subclass defs
            if isinstance(node, ParameterNode):
                # load parameter value from intermediate file
                values[node.get_name()] = open(node.input_node.filename,"rU").read().strip()
            else:
                if isinstance(node.input_node, FileNode):
                    values[node.get_name()] = '"'+node.input_node.filename+'"'
                elif isinstance(node.input_node, FolderNode):
                    values[node.get_name()] = '"'+node.input_node.foldername+'"'
                elif isinstance(node.input_node, SharedInputNode):
                    try:
                        values[node.get_name()] = '"'+node.input_node.input_node.foldername+'"'
                    except AttributeError:
                        values[node.get_name()] = '"'+node.input_node.input_node.filename+'"'


            for outnode in node.output_nodes:
                if isinstance(outnode, FileNode):
                    values[node.get_name()] = '"'+outnode.filename+'"'
                if isinstance(outnode, FolderNode):
                    values[node.get_name()] = '"'+outnode.foldername+'"'

        # allow simple named parameters        
        for k in self.__dict__:
            o = self.__dict__[k]
            if not isinstance(o, Node):
                values[k] = '"{}"'.format(o)

        try:
            #print("command: ", command)
            #print("values: ", values)
            formatted = command.format(**values)
        except KeyError as err:
            msg = "Error: for component {}, {} node not linked to a filename or parameter,"
            msg += " or that node name does not exist."
            msg = msg.format(self.get_name(),
                             err)
            raise KeyError(msg)

        return formatted

    def collect_component_nodes(self,
                                name=None,
                                parent_name=None):
        """
        Return a non-redundant list of all ComponentNodes contained
        in this component. If a name is provided, match it, allowing
        for wildcard matching on the ends.

        """
        internal = []
        for c in self.internal_components:
            internal += c.collect_component_nodes(name=name)
        nodes = self.get_component_nodes()
        if name is not None:
            # allow simple wildcard matching
            if '*' not in name:
                nodes = [n for n in nodes if n.get_name()==name]
            elif name.startswith('*') and name.endswith('*'):
                nodes = [n for n in nodes if name[1:-1] in n.get_name()]
            elif name.startswith('*'):
                nodes = [n for n in nodes if n.get_name().endswith(name[1:])]
            elif name.endswith('*'):
                nodes = [n for n in nodes if n.get_name().startswith(name[:-1])]
            else:
                raise RuntimeError("Error: collect_component_nodes() does not support matching names with internal wildcards.")
        if parent_name is not None:
            if '*' not in parent_name:
                nodes = [n for n in nodes if n.parent_component.get_name()==name]
            elif parent_name.startswith('*') and parent_name.endswith('*'):
                nodes = [n for n in nodes if parent_name[1:-1] in n.parent_component.get_name()]
            elif parent_name.startswith('*'):
                nodes = [n for n in nodes if n.parent_component.get_name().endswith(parent_name[1:])]
            elif parent_name.endswith('*'):
                nodes = [n for n in nodes if n.parent_component.get_name().startswith(parent_name[:-1])]
            else:
                raise RuntimeError("Error: collect_component_nodes() does not support matching names with internal wildcards.")

        return list(set(nodes + internal))

    def collect_components(self):
        """
        Return a list of this component + all components contained in this one.

        """
        internal = []
        for c in self.internal_components:
            internal += c.collect_components()
        return [self] + internal

    # FIXME: support wildcard matching mid-level components by name

    def collect_low_level_components(self,
                                     name=None):
        """
        Return a list of all low level components contained in this one
        (i.e. those that actually execute a subprocess). If a name is
        provided, match it, allowing for wildcards on either end.

        """
        if len(self.internal_components) == 0:
            # print(str(self))
            selfname = self.get_name()
            if name is not None:
                # allow simple wildcard matching
                if '*' not in name:
                    if selfname==name:
                        return [self]
                    else:
                        return []
                elif name.startswith('*') and name.endswith('*'):
                    if name[1:-1] in selfname:
                        return [self]
                    else:
                        return []
                elif name.startswith('*'):
                    if selfname.endswith(name[1:]):
                        return [self]
                    else:
                        return []
                elif name.endswith('*'):
                    if selfname.startswith(name[:-1]):
                        return [self]
                    else:
                        return []
                else:
                    raise RuntimeError("Error: collect_low_level_components() does not support matching names with internal wildcards.")
            else:
                return [self]
        else:
            internal = []
            for c in self.internal_components:
                internal.extend(c.collect_low_level_components(name=name))
            return internal

    def collect_edges(self):
        edges = []
        nodes = self.collect_component_nodes()
        for from_node in nodes:
            for to_node in from_node.output_nodes:
                edge = Edge(from_node,
                            to_node)
                if edge not in edges:
                    edges.append(edge)
            if from_node.input_node is not None:
                edge = Edge(from_node.input_node,
                            from_node)
                if edge not in edges:
                    edges.append(edge)
        return edges

    def __str__(self):
        s = type(self).__name__
        if self.name is not None:
            s += ' "{}"'.format(self.name)
        else:
            s += ' "{}"'.format(self.__class__.__name__)
        s += ' {}'.format(self.id)
        return s

    def read_stdout(self):
        try:
            f = open(self.stdout.output_nodes[0].filename, 'rU')
            return f.read()
        except AttributeError:
            return ""

    def read_stderr(self):
        try:
            f = open(self.stderr.output_nodes[0].filename, 'rU')
            return f.read()
        except AttributeError:
            return ""

    def cmd(self):
        """
        Override in derived classes. Returns a string to be run
        through a commandline shell.

        Input or output filepaths can be referred to by node
        name within curly braces. For example, if this component
        had an input node with name "input" and an output node
        with name "output":
            `tail {input} > {output}`

        Parameters may also be referenced in the same way.
        """
        cmd = 'echo "No command specified for component {}.'
        cmd += ' Override cmd() in derived Component class."'
        cmd += ' 1>&2'
        cmd = cmd.format(self.get_name())
        return cmd

    def proc_status(self):
        """
        Check proc returncode. Most derived classes will use this
        default, but some may need special handling. Broken pipe
        errors are ignored, since those could be caused by
        downstream module failures.

        Returns:
         "running" - process not yet terminated
         "success" - process terminated with return code 0
         "failed" - process terminated with error (excluding broken pipes)
         "terminated" - process terminated with any error (including broken pipes)
        """

        # only call this func on components with an associated subprocess
        assert self.proc is not None

        poll = self.proc.poll() # if this returns None, process has not terminated
        r = self.proc.returncode # it seems possible to sometimes get a misleading 0 returncode and 
                                 # non-None poll before proc fully terminates. Not sure how this is
                                 # possible or where the returncode change is occuring 

        # return code -11 usually indicates SEGFAULT
        # STAR returns 139 on SEGFAULT
        # 144, 141 indicate broken pipe errors
        # negative codes usually indicate termination signal
        if r == -11:
            # write seg fault message to proc stderr file (in case proc isn't a shell
            # wrapper that already handles this)
            try:
                f = open(self.stderr.output_nodes[0].filename, 'w')
                f.write("\nSegmentation fault\n")
                f.close()
            except AttributeError:
                pass
        status = "terminated"
        if r is None:
            status = "running"
        elif r == 0:
            status = "success"
        elif r != 144 and r != 141 and (r > 0 or r == -11):
            status = "failed"
        else:
            # "terminated"
            pass
        return status

    def format_stream(self,
                      stream="stdout",
                      tail=500):
        term_width, term_height = shutil.get_terminal_size()
        s = " {} ".format(stream)
        '╴'
        '╶'
        '╸'
        '╺'
        '·'
        '░'
        #s += '░' * (term_width - len(s)) + '\n'
        s += '\n'
        s += '╴' * term_width + '\n'
        try:
            if stream == "stdout":
                s += '\n'.join(self.read_stdout().splitlines()[-tail:]) + "\n"
            elif stream == "stderr":
                s += '\n'.join(self.read_stderr().splitlines()[-tail:]) + "\n"
        except FileNotFoundError:
            s += "\n"
        return s

    def format_stdout(self, tail=500):
        return self.format_stream(stream="stdout", tail=tail)

    def format_stderr(self, tail=500):
        return self.format_stream(stream="stderr", tail=tail)

    def format_proc_status(self):
        s = "\nComponent \"" + self.get_name() + "\""
        if self.assoc_rna is not None or self.assoc_sample is not None:
            s += " ("
            if self.assoc_rna is not None:
                s += "RNA:" + self.assoc_rna
            if self.assoc_sample is not None:
                if self.assoc_rna is not None:
                    s += ' '
                s += "sample:" + self.assoc_sample
            s += ')'
        if self.proc is None:
            s += " no associated process\n"
        else:
            s += " status: " + self.proc_status() + " (return code {})".format(self.proc.returncode) + "\n"
        return s

    def get_error(self,
                  verbose=False):
        """

        Returns:
            str: subprocess error message, formatted with pipeline
                 component name indicated
        """
        assert self.proc is not None
        #assert self.proc_status() in ["failed", "terminated"]

        term_width, term_height = shutil.get_terminal_size()
        s = "\nERROR: Component \"" + self.get_name() + "\""
        if verbose:
            s += " id:{}".format(self.id)
        if self.assoc_rna is not None or self.assoc_sample is not None:
            s += " ("
            if self.assoc_rna is not None:
                s += "RNA:"+self.assoc_rna
            if self.assoc_sample is not None:
                if self.assoc_rna is not None:
                    s += ' '
                s += "sample:"+self.assoc_sample
            s += ')'
        s += " failed, giving the following error message:"
        s += '=' * (term_width - len(s))
        if "Bowtie" not in self.get_name():
            s += "\n" + self.read_stderr()
        else:
            # bowtie2 wrapper redirects stderr into stdout, so read that instead
            # FIXME: make a component member function that handles this or figure out 
            #        a better way to wrap bowtie2 without redirecting stderr
            s += "\n" + self.read_stdout()
        s += '=' * term_width + "\n"
        return s

    def after_run_message(self):
        """
        Override to display a message in log and to the user
        after a particular component is run.
        Example: `return self.read_stdout()`
        """
        return None

    def start_process(self,
                      verbose=False):
        # only call this func for lowest-level components
        assert len(self.internal_components) == 0
        if False:
            print("About to start process for {}...".format(self.get_name()))
        kwargs = {"shell":True,
                  "executable":"/bin/bash",
                  "preexec_fn": os.setsid}
        try:
            kwargs["stdout"] = open(self.stdout.output_nodes[0].filename, "w")
        except AttributeError:
            pass
        try:
            kwargs["stderr"] = open(self.stderr.output_nodes[0].filename, "w")
        except AttributeError:
            pass

        formatted_cmd = self.format_command(self.cmd())
        if verbose:
            print('\n\n'+formatted_cmd+'\n\nfrom inside dir\n\n{}\n\n'.format(os.getcwd()))

        self.proc = sp.Popen(formatted_cmd,
                             **kwargs)

        if False:
            print("...successfully started")


    def get_connected_output_components(self):
        """
        Get immediately connected downstream Components, even
        if a FileNode is in between.

        """
        connected_components = []
        for outnode in self.output_nodes:
            for connected_node in outnode.output_nodes:
                if isinstance(connected_node, ComponentNode):
                    connected_components.append(connected_node.parent_component)
                else:
                    for second_node in connected_node.output_nodes:
                        if isinstance(second_node, ComponentNode):
                            connected_components.append(second_node.parent_component)

        connected_components = list(set(connected_components))
        return connected_components

    def terminate(self):
        """
        Send SIGTERM to running subprocesses.

        """
        if self.proc is not None and self.proc.poll() is None:
            pid = self.proc.pid
            try:
                os.killpg(os.getpgid(pid), signal.SIGTERM)
            except ProcessLookupError:
                pass

    def kill(self):
        """
        Send SIGKILL to running subprocesses.

        """
        if self.proc is not None and self.proc.poll() is None:
            pid = self.proc.pid
            try:
                os.killpg(os.getpgid(pid), signal.SIGKILL)
            except ProcessLookupError:
                pass


