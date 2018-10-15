# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
import sys
import shutil
from time import sleep
from timeit import default_timer as timer
import traceback
from collections import deque
from pprint import pprint

from pyshapemap.component import *
from pyshapemap.components import Mangler, ProgressMonitor
from pyshapemap.connect import connect, disconnect
from pyshapemap.flowchart import draw_flowchart
from pyshapemap.util import timestamp, format_message, non_block_read

class Pipeline(Component):
    
    #TODO: add a pipeline validator: check connections are two-sided, check no loops, check no loops in parent structure, etc.

    def __init__(self, **kwargs):
        self.verbose = False
        super().__init__(**kwargs)


    # - Ideally Want same filenames and locations for split-to-disk files as
    #   intermediate files when run in serial mode (not critical)

    def collect_entry_components(self,
                                 skip_run_order_flagged=False):
        # FIXME: unit test
        # FIXME: might be causing problems with SharedInputNode
        comps = self.collect_low_level_components()
        # entry components have no components connected to their input nodes
        entry_comps = []
        for c in comps:
            # skip components whose run order has already been set (used by set_run_order())
            if skip_run_order_flagged:
                try:
                    if c.run_order is not None:
                        continue
                except AttributeError:
                    pass
            found_connected = False
            for node in c.input_nodes:
                if node.input_node is not None:
                    connected_comp = None
                    if isinstance(node.input_node,
                                  OutputNode):
                        # directly connected to an output node of another Component
                        connected_comp = node.input_node.parent_component
                    elif isinstance(node.input_node,
                                    (FileNode, SharedInputNode)):
                        # connected to a FileNode, possibly a PipeNode
                        # - Check if there is a connected ComponentNode one step further
                        if node.input_node.input_node is not None:
                            connected_comp = node.input_node.input_node.parent_component
                    if connected_comp is not None:
                        if skip_run_order_flagged:
                            try:
                                if connected_comp.run_order is None:
                                    found_connected = True
                            except AttributeError:
                                found_connected = True
                        else:
                            found_connected = True
                if found_connected:
                    break
            if not found_connected:
                entry_comps.append(c)
        return entry_comps

    def sort_components(self,
                        components):
        """
        Sort components by location in pipeline hierarchy

        """
        # First tag every pipeline component with a pipeline location
        def recurse(o, tree_path=None):
            if tree_path is None:
                tree_path = []
            assert isinstance(o, Component)
            #print("Component {} ({}) at {}".format(o.name,
            #                                    o.id,
            #                                    tree_path))
            o.pipeline_location = tree_path
            if len(o.internal_components) == 0:
                return
            for i in range(len(o.internal_components)):
                m = o.internal_components[i]
                updated_tree_path = tree_path + [i]
                recurse(m, tree_path=updated_tree_path)
        recurse(self)

        # FIXME: add unit tests for this func
        sorted_components = sorted(components, key=lambda x: x.pipeline_location)
        return sorted_components

    def collect_parallel_components(self,
                                    component):
        parallel_components = [component]
        # traverse pipeline connections to find components that can be
        # run in parallel with the one given
        # FIXME: this does not correctly handle components with mixed parallel/serial outputs
        touched = []
        def traverse(o):
            nonlocal parallel_components
            nonlocal touched
            if o is None:
                return
            if o not in touched:
                touched.append(o)
            else:
                return
            if isinstance(o, Component):
                if o not in parallel_components:
                    parallel_components.append(o)
                for node in o.get_component_nodes():
                    traverse(node)
            if isinstance(o, ComponentNode):
                if not o.parallel:
                    return
                traverse(o.input_node) 
                for node in o.output_nodes:
                    traverse(node)
                traverse(o.parent_component)
            elif isinstance(o, PipeNode):
                traverse(o.input_node)
                for node in o.output_nodes:
                    traverse(node)
            else: # FileNode
                return
        traverse(component)
        return parallel_components

    def map_pipeline_tree(self,
                          get_return_codes=False):
        """
        Recurse all components in pipeline hierarchy,
        and keep track of tree walk for each low-level component
        (i.e. identify the location of every endpoint)

        Args:
            pipeline:

        Returns:
            List of low-level Component ids and list of
            traversal paths to Components

        """

        names = []
        ids = []
        paths = []
        return_codes = []

        def traverse(o, tree_path=""):
            assert isinstance(o, Component)
            if len(o.internal_components)==0:
                tree_path = tree_path[:-1] # strip trailing '.'
                #sys.stdout.write("Module {} ({}) at {}\n".format(o.name,
                #                                                 o.id,
                #                                                 tree_path))
                names.append(o.get_name())
                ids.append(o.id)
                paths.append(tree_path)
                if get_return_codes:
                    if ( len(o.internal_components) == 0
                         and o.proc is not None ):
                        return_codes.append(o.proc.returncode)
                    else:
                        return_codes.append("")
                    return
            for i in range(len(o.internal_components)):
                m = o.internal_components[i]
                updated_tree_path = tree_path + "{}.".format(i)
                traverse(m, tree_path = updated_tree_path)

        traverse(self)

        if get_return_codes:
            return names, ids, paths, return_codes
        return names, ids, paths

    def get_component_at_location(self,
                                  location):
        """
        Get a specific Component within the pipeline hierarchy
        given a location string.

        Args:
            pipeline:
            location:  '.' separated string of zero-based int indices.
                        For example: "3.0" indicates the 1st component in
                        the 4th component in the main pipeline

        """
        assert isinstance(location, str)
        indices = [int(x) for x in location.split('.')]
        o = self
        for i in indices:
            o = o.internal_components[i]
        return o

    def calc_run_order(self,
                       run_order=1,
                       serial_mode=False):
        """
        Analyze the pipeline graph to determine what low-level Components
        can be run in parallel, and in what order to run processes overall.
        Set the "run_order" property for each low-level Component.

        """

        # reset component run_orders (allows this to be run more than once
        # on the same pipeline)
        if run_order==1:
            for c in self.collect_low_level_components():
                c.run_order = None

        comps = self.collect_entry_components(skip_run_order_flagged=True)

        if comps is None or len(comps) == 0:
            return
        # sort by location in pipeline hierarchy
        comps = self.sort_components(comps)
        starting_comp = comps[0]
        if serial_mode:
            starting_comp.run_order = run_order
        else:
            # collect components that can be run in parallel with
            # this one
            parallel_components = self.collect_parallel_components(starting_comp)
            #print("Setting run_order to {} for components:".format(run_order))
            #for c in sorted(parallel_components, key=lambda x: x.pipeline_location):
            #    print(" {}".format(c.get_name()))
            for c in parallel_components:
                c.run_order = run_order
        updated_run_order = run_order + 1
        # continue until all low-level components have a run_order
        self.calc_run_order(run_order=updated_run_order,
                            serial_mode=serial_mode)

    def get_run_group(self,
                      run_order):
        """
        Get a list of all low-level components with the given
        run order group.

        """
        assert isinstance(run_order, int)
        run_group = []
        for comp in self.collect_low_level_components():
            try:
                ro = comp.run_order
            except AttributeError:
                raise RuntimeError("ERROR: must calc_run_order() before attempting to run pipeline")
            if ro == run_order:
                run_group.append(comp)
        run_group = self.sort_components(run_group)
        return run_group

    def find_enclosing_component(self,
                                 components):
        """
        Find a component (if any) that contains all the
        provided components and no others.

        """
        # TODO: unit test
        # FIXME: doesn't seem to be working, at least for finding Sample instances
        comp_set = set(components)
        for component in self.collect_components():
            contained_components = component.collect_low_level_components()
            if set(contained_components) == comp_set:
                return component
        return None

    def add_output_nodes(self):
        """
        Add any missing output FileNodes and FolderNodes

        """
        comps = self.collect_low_level_components()
        for c in comps:
            for node in c.output_nodes:
                if len(node.output_nodes) == 0:
                    if node.isfolder:
                        f_node = FolderNode(error_on_existing=node.error_on_existing,
                                            make_parent=node.make_parent)
                    else:
                        f_node = FileNode()
                    # TODO: should FileNodes have a parent_component set?
                    connect(node, f_node)


    def add_intermediate_nodes(self):
        """
        Add intermediate FileNodes, FolderNodes, or PipeNodes at connections between ComponentNodes
        flagged with the same execution group.

        """

        err_msg = "ERROR: Component run order not set. Must calc_run_order() "
        err_msg += "before add_intermediate_filenodes()"

        comps = self.collect_low_level_components()
        for c in comps:
            for from_node in c.output_nodes:
                assert isinstance(from_node, ComponentNode)
                output_nodes = [n for n in from_node.output_nodes if isinstance(n, ComponentNode)]
                for to_node in output_nodes: # avoid iterating over changing list
                    try:
                        ro1 = from_node.parent_component.run_order
                        ro2 = to_node.parent_component.run_order

                    except AttributeError:
                        # no run_order set for StarAlignerMixedInput, since not a low level component.
                        # - take min of run orders from parent components of "viewing" nodes
                        # FIXME: this is hacky, and doesn't support nesting
                        if isinstance(to_node, SharedInputNode):
                            ros = []
                            for onode in to_node.output_nodes:
                                ros.append(onode.parent_component.run_order)
                            ro2 = min(ros)
                        else:
                            continue
                            # don't stop execution. easier to diagnose problem
                            # if flowchart actually gets rendered.
                            #raise RuntimeError(err_msg)
                    if ro1 is None or ro2 is None:
                        continue
                        #raise RuntimeError(err_msg)

                    # use existing file node if present (so one output
                    # can be connected to multiple inputs)
                    existing_linking_node = None
                    for n in from_node.output_nodes:
                        if isinstance(n, (FileNode, FolderNode)):
                            existing_linking_node = n
                            break

                    if existing_linking_node is not None:
                        linking_node = existing_linking_node
                    elif ro1 == ro2 :
                        # components can be run in parallel, so connect with an in-memory pipe
                        linking_node = PipeNode()
                        connect(from_node, linking_node)
                    else:
                        if from_node.isfolder:
                            linking_node = FolderNode(error_on_existing=from_node.error_on_existing,
                                                      make_parent=from_node.make_parent)
                        else:
                            linking_node = FileNode()
                        connect(from_node, linking_node)

                    # remove previous connection between ComponentNodes (ugly)
                    from_node.output_nodes.pop(from_node.output_nodes.index(to_node))
                    to_node.input_node = None
                    # connect file node to component node input
                    connect(linking_node, to_node)

    def get_output_filenodes(self,
                             components=None):
        if components is None:
            node_list = self.collect_component_nodes()
        else:
            node_list = []
            for c in components:
                for node in c.output_nodes:
                    if isinstance(node, ComponentNode):
                        if node not in node_list:
                            node_list.append(node)
        filenodes = []
        for comp_node in node_list:
            for node in comp_node.output_nodes:
                if isinstance(node, FileNode):
                    if node not in filenodes:
                        filenodes.append(node)
        return filenodes

    def get_output_foldernodes(self,
                             components=None):
        if components is None:
            node_list = self.collect_component_nodes()
        else:
            node_list = []
            for c in components:
                for node in c.output_nodes:
                    if isinstance(node, ComponentNode):
                        if node not in node_list:
                            node_list.append(node)
        foldernodes = []
        for comp_node in node_list:
            for node in comp_node.output_nodes:
                if isinstance(node, FolderNode):
                    if node not in foldernodes:
                        foldernodes.append(node)
        return foldernodes

    def gen_filenames(self,
                      path=""):
        """
        Set any currently blank FileNode filenames and FolderNode foldernames
        to sensible names based on given folder, pipeline hierarchy and node name

        """
        def gen_name(node):
            extension = node.get_extension()
            parent = node.input_node.parent_component
            parents = []
            depth_limit = 500
            i = 0
            print_parent = False
            msg = ""
            while parent is not None:
                parents = [parent.get_name()]+parents
                parent = parent.parent_component
                if print_parent:
                    msg += "\n\"{}\" parent is".format(parent.get_name())
                i += 1
                if i > depth_limit:
                    print_parent = True
                if i > depth_limit + 5:
                    raise RuntimeError("Error: apparent loop in parent_component(s) for connected node \"{}\" with immediate parent \"{}\"\n".format(node.input_node.get_name(),
                                                                                                                                                     node.input_node.parent_component.get_name())+msg+". . .")
            name = os.path.join(path,
                                *parents)
            name = os.path.join(name,
                                '_'.join(parents) + '_' +
                                node.input_node.get_name())
            if len(extension)>0:
                name += "." + extension
            return name

        for node in self.get_output_filenodes():
            if node.filename is not None:
                continue
            node.filename = gen_name(node)

        for node in self.get_output_foldernodes():
            if node.foldername is not None:
                continue
            node.foldername = gen_name(node)

    def make_paths(self,
                   components=None):
        """
        Create enclosing paths to output and intermediate files if needed
        If a list of components is provided, only make paths for files
        associated with those components.
        """
        for node in list(set(self.get_output_filenodes(components=components) + \
                    self.get_output_foldernodes(components=components))):
            node.make_path()


    def remove_existing_pipes(self,
                              components=None):
        """
        Remove named pipes if present (otherwise there can be problems
        with repeated runs in the same folder if a process attempts to
        open a previously created pipe as a file).

        If a list of components is provided, only examine paths
        associated with those components.
        """
        for node in self.get_output_filenodes(components=components):
            node.remove_existing_pipe()

    def remove_existing_files(self,
                              components=None):
        for node in self.get_output_filenodes(components=components):
            node.remove_existing_file()

    def remove_existing_folders(self,
                                components=None):
        for node in self.get_output_foldernodes(components=components):
            node.remove_existing_folder()


    def make_pipes(self,
                   components=None):
        """
        Create named pipes (FIFOs) for PipeNodes
        If a list of components is provided, only make paths for files
        associated with those components.
        """
        for node in list(set(self.get_output_filenodes(components=components))):
            if isinstance(node, PipeNode):
                node.make_pipe()

    def setup(self,
              temp=None,
              serial_mode=None,
              **kwargs):
        """
        Perform actions needed to prepare pipeline
        for execution.

        """
        if temp is None:
            temp = self.temp
        if serial_mode is None:
            serial_mode = self.serial
        self.calc_run_order(serial_mode=serial_mode)
        self.add_output_nodes()
        self.add_intermediate_nodes()
        self.gen_filenames(path=temp)

    def get_terminal_components(self,
                                components=None):
        """
        Get terminal low-level pipeline components. Limit
        to given list of components if specified. Assumes
        presence of intermediate FileNodes in pipeline.

        """
        terminal_components = []
        if components is None:
            components = self.collect_low_level_components()
        # terminal components connected output nodes are connected
        # to no components in this group
        for component in components:
            found_connection = False
            for node in component.output_nodes:
                for filenode in node.output_nodes:
                    for compnode in filenode.output_nodes:
                        if compnode.parent_component in components:
                            found_connection = True
            if not found_connection:
                terminal_components.append(component)
        return terminal_components


    def get_upstream_failed_components(self,
                                       components):
        """
        Find the most upstream components in a run group
        responsible for failure.

        """
        if len(components) == 1:
            return components

        failed_components = []
        for c in components:
            if c.proc_status() == "failed":
                failed_components.append(c)

        def check_connected(from_component,
                            to_component):
            to_crawl = deque([from_component])
            visited_components = set()
            while to_crawl:
                current = to_crawl.popleft()
                if current in visited_components:
                    continue
                if current is to_component:
                    return True
                visited_components.add(current)
                component_children = []
                for connected_comp in current.get_connected_output_components():
                    component_children.append(connected_comp)
                component_children = set(component_children)
                to_crawl.extend(component_children - visited_components)
            return False

        # if any downstream traversal connects two failed components,
        # don't report the downstream component
        keep_component = [True for c in failed_components]
        for i in range(len(failed_components)):
            for j in range(len(failed_components)):
                if j == i:
                    continue
                if check_connected(failed_components[i],
                                   failed_components[j]):
                    keep_component[j] = False
        upstream_failed_components = []
        for n in range(len(keep_component)):
            if keep_component[n]:
                upstream_failed_components.append(failed_components[n])

        return upstream_failed_components


    def get_error_message(self,
                          components=None,
                          verbose_ids=False,
                          verbose=False,
                          tail=50):
        """
        Get error message after a running process group fails
        """
        term_width, term_height = shutil.get_terminal_size()
        if components is None:
            # if no process group specified, use the components that
            # were last run
            components = self.get_run_group(self.current_run_group_index)
        failed_components = self.get_upstream_failed_components(components)
        msg = ""
        for c in failed_components:
            msg += c.get_error(verbose=verbose_ids)
        if verbose:
            for c in components:
                msg += '▀' * term_width
                msg += c.format_proc_status()
                msg += '─' * term_width + "\n"
                msg += c.format_stdout(tail=tail)
                msg += '─' * term_width + "\n"
                msg += c.format_stderr(tail=tail)
                msg += '▄' * term_width + "\n"

        if self.hung:
            s = ' RUN EXCEEDED TIMEOUT '
            s += '◷ ' * int((term_width - len(s))/2) + '\n'
            msg += s
        return msg


    # def get_all_messages(self,
    #                      components=None,
    #                      tail=200):
    #     term_width, term_height = shutil.get_terminal_size()
    #     if components is None:
    #         components = self.get_run_group(self.current_run_group_index)
    #     msg = ""
    #     for c in components:
    #         msg += '▀' * term_width
    #         msg += c.format_proc_status()
    #         if self.hung:
    #             s = ' RUN EXCEEDED TIMEOUT '
    #             s += '◷' * (term_width - len(s)) + '\n'
    #             msg += s
    #         msg += '─' * term_width + "\n"
    #         msg += c.format_stdout(tail=tail)
    #         msg += '─' * term_width + "\n"
    #         msg += c.format_stderr(tail=tail)
    #         msg += '▄' * term_width + "\n"
    #     return msg


    def run(self,
            timeout=None,
            quiet=False,
            latency=0.1,
            success_term_pause=0.05,
            fail_term_pause=0.5):
        """
        Run process groups in order specified by
        the "run_order" component property

        """

        if self.overwrite:
            self.remove_existing_pipes()
            self.remove_existing_folders()
            self.remove_existing_files()
        # create paths to files and pipes before doing
        # anything else, since existing files will cause
        # an error message if --overwrite is not used
        self.make_paths()
        self.make_pipes()

        i = 1
        success = True
        while True:
            self.current_run_group_index = i
            run_group = self.get_run_group(i)
            if len(run_group) < 1:
                # no more components to run
                break

            # find most upstream ProgressMonitor component in process group and
            #  use that for progress bar
            progmons = [c for c in run_group if isinstance(c, ProgressMonitor)]
            if len(progmons)==1:
                progmon = progmons[0]
            elif len(progmons)==0:
                progmon = None
            else:
                progmons = self.sort_components(progmons)
                progmon = progmons[0]

            # Tell the user what is being run
            # - if there is an enclosing component that describes
            #   this run group, display that name
            # - list low-level component names
            enclosing_component = None
            if not quiet:
                if len(run_group) == 1:
                    s = "Running {} at {} . . .\n".format(run_group[0].get_name(),
                                                          timestamp())
                else:
                    enclosing_component = self.find_enclosing_component(run_group)
                    if enclosing_component is not None:
                        s = "Running {} at {} . . .\n".format(enclosing_component.get_name(),
                                                              timestamp())
                    else:
                        s = "Running process group {} at {} . . .\n".format(i,
                                                                            timestamp())
                    s += "  Including these components:\n"
                sys.stdout.write(s)

            terminal_components = self.get_terminal_components(components=run_group)


            for component in run_group:
                if not quiet:
                    if len(run_group) > 1:
                        sys.stdout.write("    {} . . .".format(component.get_name()))
                component.start_process(verbose=self.verbose)
                if not quiet:
                    if len(run_group) > 1:
                        sys.stdout.write(" started at {}\n".format(timestamp()))

            # Loop until all terminal (final) components stop
            # or any component stops with a non-zero return code
            # or timeout
            start_time = timer()
            while True:
                if timeout is not None:
                    if timer() - start_time > timeout:
                        if not quiet:
                            sys.stderr.write("Run hung or timed out.\n")
                        self.hung = True
                        for c in self.collect_low_level_components():
                            c.hung = True
                        success = False
                        break_flag = True
                        break
                sleep(latency)

                break_flag = False

                if all(c.proc_status() == "success"
                       for c in terminal_components):
                    break_flag = True
                if any(c.proc_status() in ["failed", "terminated"]
                       for c in run_group):
                    success = False
                    break_flag = True

                if break_flag:
                    # TODO: is there any way we can force stream flushes before terminating processes?
                    #        - for debugging, it's nice to have the line present in an output file that
                    #          actually caused an error in a downstream process
                    #        - could SIGTERM, then escalate to SIGKILL
                    # - for now, just waiting a bit
                    if success:
                        sleep(success_term_pause)
                    else:
                        sleep(fail_term_pause) # wait a bit longer on failure so other procs can terminate, write errors
                                               # and buffers to disk, etc.
                    break

                # Progress bar handling, if enclosing_component has
                # an associated pipeviewer component
                # FIXME: newline doesn't always print after progress bar complete
                if ( progmon is not None
                     and not quiet ):
                    o = non_block_read(progmon.proc.stderr)

                    if o is not None:
                        s = o.decode("utf-8")
                        level = 2
                        if len(s) > 1:
                            lev = ''.join([' ' for x in range(level)])
                            s = s.strip()
                            term_width, term_height = shutil.get_terminal_size()
                            pad = ' ' * (term_width - len(s) - level)
                            sys.stdout.write("\r" + lev + s + pad)
                            of = open(progmon.stderr.output_nodes[0].filename, "w")
                            of.write(s + "\n")
                            of.close()

            if not quiet:
                if success:
                    s = ""
                    # display stats/output messages for some components
                    indent = 0
                    if len(run_group) > 1:
                        indent = 2
                    for component in run_group:
                        msg = component.after_run_message()
                        title = component.get_name()
                        if component.assoc_rna:
                            title += " (RNA: {})".format(component.assoc_rna)
                        if component.assoc_sample:
                            title += " (sample: {})".format(component.assoc_sample)
                        msg = format_message(title,
                                             indent,
                                             msg)
                        s += msg
                    s += ". . . done at {}\n".format(timestamp())
                else:
                    s = self.get_error_message(run_group, verbose=self.verbose)
                sys.stdout.write(s)
            if not success:
                self.clean_up()
                break
            i += 1
        return success

    def get_failed_targets(self):
        """
        Check if quality control checks failed for any targets
        """
        qc_comps = self.collect_low_level_components(name="RenderFigures")
        failed_rnas = []
        for comp in qc_comps:
            s = comp.read_stdout()
            if "FAIL" in s:
                failed_rnas.append(comp.assoc_rna)
        return failed_rnas

    def any_star_segfault(self):
        any_fault = False
        for c in self.collect_low_level_components("StarAligner*"):
            if c.proc.returncode == 139:
                any_fault = True
                break
        return any_fault

    def terminate(self):
        for c in self.collect_low_level_components():
            c.terminate()

    def kill(self):
        for c in self.collect_low_level_components():
            c.kill()

    def clean_up(self,
                 latency=1):
        '''
        Send SIGTERM to child processes, wait 1 second, then
        send SIGKILL
        '''
        self.terminate()
        sleep(latency)
        self.kill()



# TODO: maybe move this to "test" folder
def run_fail_test(pipeline=None,
                  failing_module_loc=None,
                  failing_module_name=None,
                  quiet=True,
                  latency=0.1,
                  timeout=30,
                  success_term_pause=0.05,
                  fail_term_pause=0.5):
    p = failing_module_loc
    name = failing_module_name

    print("[ RUN      ] Fail test for component at {} ({})".format(p,
                                                                   name))

    # TODO: allow use of existing aligner index (this should dramatically speed up overall test execution)
    passed = False
    success = 1
    err = ""
    try:

        # locate component immediately containing selected endpoint component
        indices = [int(x) for x in p.split('.')]
        parent_component = pipeline
        for i in indices[:-1]: # only go down to the next-to-last element
            parent_component = parent_component.internal_components[i]

        # locate target component
        component_index = indices[-1]
        selected_component = parent_component.internal_components[component_index]

        # insert an error injector in front of each input
        input_nodes = selected_component.input_nodes
        added_count = 0
        for node in input_nodes:
            # skip nodes that input from folders, not files
            # (don't have a good solution for mangling these inputs)
            if node.get_name() == "index":
                continue

            if added_count > 0:
                name = "Mangler{}".format(added_count+1)
            else:
                name = "Mangler"
            mangler = Mangler(name=name)

            added_count += 1

            # rewire connections
            prev_node = node.input_node
            if prev_node is not None:
                disconnect(prev_node, node)
                connect(prev_node, mangler.input)
            connect(mangler.output, node)

            # update component hierarchy
            parent_component.insert(mangler, selected_component)

        if True:
            draw_flowchart(pipeline,
                           pipeline.flowchart_path[:-4]+"_pre-setup.svg",
                           path=pipeline.temp)

        # run setup() to generate missing intermediate files and/or pipes
        # Note: when running normally (not in a failure test) this is done automatically
        #       in build_pipeline()
        pipeline.setup()

        if True:
            draw_flowchart(pipeline,
                           pipeline.flowchart_path,
                           path=pipeline.temp)

        success = pipeline.run(timeout=timeout,
                               quiet=quiet,
                               latency=latency,
                               success_term_pause=success_term_pause,
                               fail_term_pause=fail_term_pause)
        # ^ fails to detect bowtie2 exit with latency 0, 0.001, 0.01, but 0.1 works

        if not success:
            err += pipeline.get_error_message(verbose_ids=True,
                                              verbose=pipeline.verbose).strip()
            fail_mod_id = None
            if len(err)>1 and err.startswith("ERROR:"):
                for s in err.splitlines()[0].split():
                    if s.startswith("id:"):
                        fail_mod_id = s[3:]
                        break
            if fail_mod_id == selected_component.id:
                passed = True
        pipeline.clean_up(latency=0) # quickly kill any hung subprocesses

    finally:
        try:
            pipeline.clean_up()
        except (NameError, AttributeError):
            pass

    return passed, success, err

