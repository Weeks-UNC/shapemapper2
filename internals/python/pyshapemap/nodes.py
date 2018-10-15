# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
import stat
import shutil

from pyshapemap.util import get_extension, rand_id

# TODO: add functions to simplify getting filenames from ComponentNodes connected node?
# FIXME: overwrite warning does not currently apply to "files" that represent a prefix passed to a process creating multiple output files

name_collision_msg = "Remove previous output folders before rerunning, or"
name_collision_msg += " use the --overwrite option to allow ShapeMapper"
name_collision_msg += " to remove existing files automatically. If running"
name_collision_msg += " multiple shapemapper instances from the same folder,"
name_collision_msg += " give each run a unique name with --name."

class Node():
    def __init__(self, **kwargs):
        kw = {"name": None,
              "input_node": None,
              "output_nodes": [],
              "id": rand_id(),
              "assoc_rna": None,
              "assoc_sample": None}
        # Nodes that are not PipeNodes
        # can have multiple nodes connected
        # in the outgoing direction
        kw.update(kwargs)
        self.__dict__.update(kw)

    def __str__(self):
        s = self.__class__.__name__
        s += ' "{}"'.format(self.name)
        s += ' {}'.format(self.id)
        return s

    def get_name(self):
        name = self.__class__.__name__
        if self.name is not None:
            name = self.name
        return name



class FileNode(Node):
    def __init__(self,
                 gv_props=None,
                 **kwargs):
        kw = {"filename": None}
        kw.update(kwargs)
        self.gv_props = {"fillcolor": "yellow",
                         "style": 'filled',
                         "shape": 'box',
                         "fontsize": '18',
                         "mono": False,
                         "margin":"0.05,0.05"}
        if gv_props is not None:
            self.gv_props.update(gv_props)
        super().__init__(**kw)

    def get_extension(self):
        # if a filename is set, use that extension,
        # otherwise attempt to use input node's extension
        if (self.filename is not None and
                    len(self.filename) > 0):
            return get_extension(self.filename)
        elif self.input_node is not None:
            return self.input_node.get_extension()
        else:
            return ""

    def make_path(self):
        """
        Create enclosing path to file if file does not exist

        """
        assert isinstance(self.filename, str)
        if not os.path.exists(self.filename):
            folder, filename = os.path.split(self.filename)
            if len(folder) > 0:
                try:
                    os.makedirs(folder, exist_ok=True)
                except NotADirectoryError:
                    raise RuntimeError("File {} already exists. ".format(folder)+name_collision_msg)
        else:
            raise RuntimeError("File {} already exists. ".format(self.filename)+name_collision_msg)

    def remove_existing_pipe(self):
        """
        If a named pipe exists at filename, delete it. (Prevents issues with
        repeated runs in the same directory)

        """
        if os.path.exists(self.filename):
            if stat.S_ISFIFO(os.stat(self.filename).st_mode):
                os.remove(self.filename)

    def remove_existing_file(self):
        if os.path.exists(self.filename):
            if stat.S_ISREG(os.stat(self.filename).st_mode):
                os.remove(self.filename)


class FolderNode(Node):
    def __init__(self,
                 gv_props=None,
                 **kwargs):
        kw = {"foldername": None,
              "extension": "",
              "error_on_existing": False, # exit with error if folder already exists (needed for STAR aligner)
              "make_parent": False} # only create folder path up to parent (needed for STAR aligner)
        kw.update(kwargs)
        self.gv_props = {"fillcolor": "yellow",
                         "style": 'filled',
                         "shape": 'box',
                         "fontsize": '18',
                         "mono": False,
                         "margin":"0.05,0.05"}
        if gv_props is not None:
            self.gv_props.update(gv_props)
        super().__init__(**kw)

    def get_extension(self):
        return ""

    def make_path(self):
        assert isinstance(self.foldername, str)
        if not os.path.isdir(self.foldername):
            # make path up to but not including final path component if option set
            if self.make_parent is not None and self.make_parent:
                folder, filename = os.path.split(self.foldername)
            else:
                folder = self.foldername
            if len(folder) > 0:
                try:
                    os.makedirs(folder, exist_ok=True)
                except NotADirectoryError:
                    raise RuntimeError("File {} already exists. ".format(folder)+name_collision_msg)
        elif self.error_on_existing is not None and self.error_on_existing:
            raise RuntimeError("Folder {} already exists. ".format(self.foldername)+name_collision_msg)

    def remove_existing_folder(self):
        if os.path.exists(self.foldername):
            if stat.S_ISDIR(os.stat(self.foldername).st_mode):
                shutil.rmtree(self.foldername, ignore_errors=True)


class PipeNode(FileNode):
    def __init__(self, **kwargs):
        gv_props = {"color": 'lightyellow2',
                    "fillcolor": 'lightyellow2',
                    "style": 'filled',
                    "shape": 'box',
                    "fontsize": '8',
                    "mono": False}
        super().__init__(gv_props=gv_props,
                         **kwargs)

    def make_pipe(self):
        assert isinstance(self.filename, str)
        try:
            os.mkfifo(self.filename)
        except FileExistsError:
            raise RuntimeError("File {} already exists. ".format(self.filename)+name_collision_msg)


class ComponentNode(Node):
    """
    A Node with an associated "parent" Component
    """

    def __init__(self,
                 gv_props=None, 
                 **kwargs):
        self.gv_props = {"color": 'gray70',
                         "style": 'filled',
                         "shape": 'box',
                         "fontsize": '14'}
        if gv_props is not None:
            self.gv_props.update(gv_props)
        kw = {"name": None,
              "parent_component": None,
              "parallel": True,
              "isfolder": False}
        kw.update(kwargs)
        if kw["isfolder"]:
            kw["parallel"] = False
        super().__init__(**kw)


class InputNode(ComponentNode):
    """
    InputNodes can have only one connected file/folder or component node

    """

    def __init__(self,
                 filename=None,
                 foldername=None,
                 **kwargs):
        kw = {"name": "input"}
        kw.update(kwargs)
        super().__init__(**kw)
        if filename is not None:
            self.set_file(filename)
        if foldername is not None:
            self.set_folder(foldername)

    def get_extension(self):
        return self.input_node.get_extension()

    def set_file(self,
                 filename):
        assert isinstance(filename, str)
        if self.input_node is not None:
            self.input_node.filename = filename
        else:
            file_node = FileNode(filename=filename)
            connect(file_node, self)

    def set_folder(self,
                   foldername):
        assert isinstance(foldername, str)
        if self.input_node is not None:
            assert isinstance(self.input_node, FolderNode)
            self.input_node.foldername = foldername
        else:
            folder_node = FolderNode(foldername=foldername)
            connect(folder_node, self)

# FIXME: it might be more consistent to have a separate OutputFolderNode class, instead of
#        the isfolder field

class OutputNode(ComponentNode):
    def __init__(self,
                 filename=None,
                 foldername=None,
                 **kwargs):
        kw = {"name": "output",
              "isfolder": False,
              "extension": "txt",
              "make_parent": False,
              "error_on_existing": False}
        kw.update(kwargs)
        super().__init__(**kw)
        if filename is not None:
            self.set_file(filename)
        if foldername is not None:
            self.set_folder(foldername)

    def get_extension(self):
        if isinstance(self.extension, str):
            if self.extension == "passthrough":
                return self.parent_component.input_nodes[0].get_extension()
            else:
                return self.extension

    # TODO: check if this function is being used anymore
    def set_file(self,
                 filename):
        assert isinstance(filename, str)
        if len(self.output_nodes)==1 and self.output_nodes[0] is not None:
            self.output_nodes[0].filename = filename
        else:
            file_node = FileNode(filename=filename)
            connect(self, file_node)

    def set_folder(self,
                   foldername):
        assert isinstance(foldername, str)
        if len(self.output_nodes)==1 and self.output_nodes[0] is not None:
            assert isinstance(self.output_nodes[0], FolderNode)
            self.output_nodes[0].foldername = foldername
        else:
            folder_node = FolderNode(foldername=foldername,
                                     error_on_existing=self.error_on_existing,
                                     make_parent=self.make_parent)
            self.output_nodes = [folder_node]


class SharedInputNode(InputNode):
    '''
    Attempting to allow two internal components to share a view of a single input node
    '''
    # - need to reexamine functions in pipeline.py such as
    #   create_intermediate_nodes, anything that traverses the workflow graph
    # - update connect() to handle linking this node to two InputNodes
    # - also update flowchart rendering
    # FIXME: filenames not getting set. check create_intermediate_nodes
    # - should probably rework to use recursion (right now probably can't nest these)

class StdinNode(InputNode):
    def __init__(self, **kwargs):
        kw = {"name": "stdin"}
        kw.update(kwargs)
        super().__init__(**kw)


class StdoutNode(OutputNode):
    def __init__(self, **kwargs):
        kw = {"name": "stdout",
              "parallel": False}
        kw.update(kwargs)
        super().__init__(**kw)


class StderrNode(OutputNode):
    def __init__(self, **kwargs):
        kw = {"name": "stderr",
              "parallel": False}
        kw.update(kwargs)
        super().__init__(**kw)


class ParameterNode(ComponentNode):
    """
    A parameter node reads a single parameter from a connected
    input file just before parent Component process execution.

    """

    # TODO: this should perhaps inherit from InputNode? 
    #       - not all InputNode methods apply though

    def __init__(self,
                 filename=None,
                 **kwargs):
        kw = {"name": "input",
              "gv_props": {"color": 'darkslategray3'}}
        kw.update(kwargs)
        super().__init__(**kw)
        if filename is not None:
            self.set_file(filename)

    def get_extension(self):
        return self.input_node.get_extension()

    def set_file(self,
                 filename):
        assert isinstance(filename, str)
        if self.input_node is not None:
            self.input_node.filename = filename
        else:
            file_node = FileNode(filename=filename)
            connect(file_node, self)


def connect(from_node, to_node):
    """
    Connect two nodes. Note that ComponentNodes or FileNodes
    that are not PipeNodes can have multiple connected nodes in the
    outgoing direction.
    """
    assert isinstance(from_node, Node)
    assert isinstance(to_node, Node)
    if (isinstance(from_node, (ComponentNode, FileNode, FolderNode)) and
        not isinstance(from_node, PipeNode)):
        from_node.output_nodes.append(to_node)
    else:
        from_node.output_nodes = [to_node]
    to_node.input_node = from_node


def connect_shared_input_nodes(from_node, to_nodes):
    """
    Connect SharedInputNode to multiple InputNodes
    """
    assert isinstance(from_node, SharedInputNode)
    assert isinstance(to_nodes, list)
    from_node.output_nodes = to_nodes
    for to_node in to_nodes:
        assert isinstance(to_node, InputNode)
        to_node.input_node = from_node


def disconnect(from_node, to_node):
    assert isinstance(from_node, Node)
    assert isinstance(to_node, Node)
    # TODO: also support reversed args for convenience?
    try:
        index = from_node.output_nodes.index(to_node)
    except ValueError:
        raise RuntimeError("Attempted to disconnect() two nodes that are not connected")
    from_node.output_nodes.pop(index)
    to_node.input_node = None


def disconnect_shared_input_nodes(from_node, to_nodes):
    """
    Disconnect SharedInputNode to multiple InputNodes
    """
    assert isinstance(from_node, SharedInputNode)
    assert isinstance(to_nodes, list)
    from_node.output_nodes = []
    for to_node in to_nodes:
        assert isinstance(to_node, InputNode)
        to_node.input_node = None


class Edge():
    def __init__(self, from_node, to_node):
        self.from_node = from_node
        self.to_node = to_node

    def __eq__(self, other):
        if (self.from_node == other.from_node and
                    self.to_node == other.to_node):
            return True
        else:
            return False

    def __str__(self):
        return "({})->({})".format(self.from_node,
                                   self.to_node)
