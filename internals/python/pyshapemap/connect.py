# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

# FIXME: figure out better way to organize these function names, so that
#        imports aren't so fragile
from pyshapemap.nodes import connect as connect_nodes
from pyshapemap.nodes import disconnect as disconnect_nodes
from pyshapemap.nodes import Node, SharedInputNode, connect_shared_input_nodes, disconnect_shared_input_nodes

from pyshapemap.component import Component


def connect(x, y):
    """
    Convenience function to link two Components with only
    a single input/output node (excluding stdout and stderr). Will also
    connect two nodes using connect() defined in nodes.py.
    """
    if isinstance(x, Node) and isinstance(y, Node):
        connect_nodes(x, y)
    elif isinstance(x, Component) and isinstance(y, Component):
        output_nodes = [node for node in x.output_nodes if node.get_name() not in ["stdout", "stderr"]]
        if len(output_nodes) != 1 or len(y.input_nodes) != 1:
            raise TypeError(
                "To connect() two Components, the first must have only one output node (excluding stdout and stderr), and the second must have only one input node.")
        connect_nodes(output_nodes[0], y.input_nodes[0])
    elif isinstance(x, SharedInputNode) and isinstance(y, list):
        connect_shared_input_nodes(x, y)
    else:
        raise TypeError("Objects to connect() must both be Nodes or Components")


def disconnect(x, y):
    if isinstance(x, Node) and isinstance(y, Node):
        disconnect_nodes(x, y)
    elif isinstance(x, Component) and isinstance(y, Component):
        output_nodes = [node for node in x.output_nodes if node.get_name() not in ["stdout", "stderr"]]
        if len(output_nodes) != 1 or len(y.input_nodes) != 1:
            raise TypeError(
                "To disconnect() two Components, the first must have only one output node (excluding stdout and stderr), and the second must have only one input node.")
        disconnect_nodes(output_nodes[0], y.input_nodes[0])
    elif isinstance(x, SharedInputNode) and isinstance(y, list):
        disconnect_shared_input_nodes(x, y)
    else:
        raise TypeError("Objects to connect() must both be Nodes or Components")
