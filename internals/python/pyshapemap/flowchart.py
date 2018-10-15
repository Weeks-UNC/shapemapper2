# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
import subprocess as sp

import pyshapemap.components as cmp


def display_image(filename):
    cmd = ["gnome-open", filename]
    sp.Popen(cmd, preexec_fn=os.setsid, stderr=sp.PIPE, stdout=sp.PIPE)


def tab(s):
    return '\n'.join(['\t' + line for line in s.splitlines()]) + '\n'

def format_path(s,
                path=None,
                replace_long_filename=False):
    """
    Break up nested folder strings vertically to prevent large horizontal boxes.
    Also replace shared starting folder with [OUT] and replace leading part of long
    filename with ". . ."
    """
    # TODO: could possibly simplify by moving representation logic to FileNode
    replaced_path_name = False
    if path is not None:
        if s.startswith(path):
            replaced_path_name = True
            s = s[len(path)+1:]
    if replace_long_filename:
        head, tail = os.path.split(s)
        name_prefix = head.replace('/','_')
        if '/' in head and len(name_prefix) > 0:
            s = s.replace(name_prefix, '. . . ')
    if replaced_path_name:
        s = "[DIR]/" + s
    return "/\\n".join(s.split('/'))

def body_str(legend, subgraphs, edges):
    s = 'strict digraph {{\n\tsplines=false;outputorder=nodesfirst;rankdir=TB;\n{}\n{}\n{}\n}}'
    s = s.format(tab(legend),
                 tab(subgraphs),
                 tab(edges))
    return s

def legend_str(path=""):
    # can use rank=same within graph to evenly space legend horizontally, but
    # can't seem to get it to appear above the main flowchart
    legend="""subgraph cluster_legend {{
margin=20;
graph [color=gray,fontname="Arial",fontsize=26,label="\\nKey",labeljust=l,shape=box,];
{{
b0 [color="white",fillcolor="white",style=filled,fontname="Arial",fontsize=14,label="ShapeMapper run from: {}",shape=box,width=2];
b01 [color="white",fillcolor="white",style=filled,fontname="Arial",fontsize=14,label="[DIR] = {}",shape=box,width=2];
b1 [color=black,fillcolor=white,fontname="Arial",fontsize=14,label=component,shape=box,style=filled,width=2];
b2 [color="gray70",fontname="Arial",fontsize=14,label="input or output node",shape=box,style=filled,width=2];
b22 [color="darkslategray3",fontname="Arial",fontsize=14,label="input parameter node",shape=box,style=filled,width=2];
b3 [fillcolor="mediumseagreen",style=filled,fontname="Arial",fontsize=14,label="input file/folder",shape=box,style=filled,width=2];
b4 [fillcolor="maroon1",style=filled,fontname="Arial",fontsize=14,label="output file/folder",shape=box,width=2];
b41 [fillcolor="yellow",style=filled,fontname="Arial",fontsize=14,label="main output file",shape=box,width=2];
b5 [fillcolor="slategray1",style=filled,fontname="Arial",fontsize=14,label="stdout or stderr file",shape=box,width=2];
b6 [color="lightyellow2",fillcolor="lightyellow2",style=filled,fontname="Arial",fontsize=14,label="pipe",shape=box,width=2];
}}
edge[style=invis,weight=1000];
b0->b01->b1->b2->b22->b3->b4->b41->b5->b6;
}}"""
    legend = legend.format(format_path(os.getcwd()),
                           path)
    return legend

def subgraph_str(internal='',
                 name='',
                 label='',
                 style='',
                 color='',
                 penwidth='',
                 shape='',
                 fontsize='',
                 margin='',
                 nodesep='',
                 ranksep=''):
    s = 'subgraph {} {{\nhandle [style=invis];{}\n}}\n'
    k = node_str('graph',
                 label=label,
                 color=color,
                 penwidth=penwidth,
                 style=style,
                 shape=shape,
                 fontsize=fontsize,
                 margin=margin,
                 nodesep=nodesep,
                 ranksep=ranksep)
    s = s.format(name, tab(k + internal))
    return s


def node_str(ID,
             label='',
             color='',
             fillcolor='',
             penwidth='',
             style='',
             shape='',
             fontsize='',
             mono=False,
             margin='',
             nodesep='',
             ranksep=''):
    # TODO: could replace most of this func with a loop over kwargs
    s = '{} [color={},'
    s += 'fillcolor={},'
    s += 'penwidth={},'
    if mono:
        s += 'fontname="Mono",'
    else:
        s += 'fontname="Arial",'  # Arial confuses l and I
    s += 'fontsize={},'
    s += 'label={},'
    s += 'shape={},'
    s += 'style={},'
    s += 'margin={},'
    s += 'nodesep={},'
    s += 'ranksep={},'
    s += '];\n'

    if not ID.isalpha():
        ID = '"{}"'.format(ID)
    if not color.isalpha():
        color = '"{}"'.format(color)
    if not fillcolor.isalpha():
        fillcolor = '"{}"'.format(fillcolor)
    if not penwidth.isalpha():
        penwidth = '"{}"'.format(penwidth)
    if not label.isalpha():
        label = '"{}"'.format(label)
    if not style.isalpha():
        style = '"{}"'.format(style)
    if not margin.isalpha():
        margin = '"{}"'.format(margin)
    if not nodesep.isalpha():
        nodesep = '"{}"'.format(nodesep)
    if not ranksep.isalpha():
        ranksep = '"{}"'.format(ranksep)
    s = s.format(ID, color, fillcolor, penwidth, fontsize, label, shape, style, margin, nodesep, ranksep)
    return s


def edge_str(ID1, ID2, style="", weight="", minlen="", color=""):
    s = '{} -> {} {};\n'
    k = ''
    if not color.isalpha():
        color = '"{}"'.format(color)
    if not style.isalpha():
        style = '"{}"'.format(style)
    if not weight.isalpha():
        weight = '"{}"'.format(weight)
    if not minlen.isalpha():
        minlen = '"{}"'.format(minlen)
    k = '[color={}, style={}, weight={}, minlen={}]'.format(color, style, weight, minlen)
    if not ID1.isalpha():
        ID1 = '"{}"'.format(ID1)
    if not ID2.isalpha():
        ID2 = '"{}"'.format(ID2)
    s = s.format(ID1, ID2, k)
    return s


def pipeline_to_dot(pipeline,
                    path="",
                    highlight_dir=None):

    assert isinstance(pipeline, cmp.Component)

    touched = []

    def traverse(o):
        nonlocal touched

        if o is None:
            return "", ""

        if o not in touched:
            touched.append(o)
        else:
            return "", ""

        if isinstance(o, cmp.Component):
            child_scope = ""
            same_scope = ""
            if len(o.internal_components) == 0:
                # attempt to preserve horizontal order of nodes
                # (will stay horizontal, but may flip L->R or R->L)
                for nodes in [o.input_nodes, o.output_nodes]:
                    for i in range(len(nodes)-1):
                        n1 = nodes[i]
                        n2 = nodes[i+1]
                        s = ""
                        s += "{\nrank=same;\n"
                        s += edge_str(n1.id,
                                      n2.id,
                                      style="invis",
                                      weight="1")
                        s += "}\n"
                        child_scope += s
                # traverse connected nodes
                for node in o.input_nodes:
                    s, p = traverse(node)
                    child_scope += s
                    same_scope += p
                for node in o.output_nodes:
                    s, p = traverse(node)
                    child_scope += s
                    same_scope += p

            # else:
            for m in o.internal_components:
                s, p = traverse(m)
                child_scope += s
                same_scope += p
            label = o.get_name()
            try:
                if o.assoc_sample is not None:
                    label += "\\n({})".format(o.assoc_sample)
                if o.assoc_rna is not None:
                    label += "\\n({})".format(o.assoc_rna)
                if o.run_order is not None:
                    label += "\\n(run group {})".format(o.run_order)
            except AttributeError:
                pass
            if o.parent_component is None:
                # TODO: add some additional run info to outermost pipeline
                pass
            # reduce the space between internal nodes and box edges
            # if this is a low-level component
            if len(o.internal_components) == 0:
                o.gv_props["margin"] = "0"
                o.gv_props["nodesep"] = "0"
                o.gv_props["ranksep"] = "0"
            s = subgraph_str(internal=child_scope,
                             name="cluster_" + o.id,
                             label=label,
                             **o.gv_props) + same_scope

            return s, ""

        elif isinstance(o, cmp.PipeNode):
            s = node_str(o.id,
                         label=format_path(o.filename,
                                           path=path,
                                           replace_long_filename=True),
                         **o.gv_props)
            return s, ""
        elif isinstance(o, (cmp.FileNode, cmp.FolderNode)):
            if o.input_node is None:
                o.gv_props["fillcolor"] = "mediumseagreen"
            else:
                o.gv_props["fillcolor"] = "maroon1"
                # color files in main output folder differently
                if highlight_dir is not None and highlight_dir != path:
                    if isinstance(o, cmp.FolderNode):
                        fname = o.foldername
                    elif isinstance(o, cmp.FileNode):
                        fname = o.filename
                    if fname is not None and fname.startswith(highlight_dir):
                        o.gv_props["fillcolor"] = "yellow"
                if o.input_node.get_name() in ["stdout", "stderr"]:
                    o.gv_props["fillcolor"] = "slategray1"
                    o.gv_props["fontsize"] = "8"
            if isinstance(o, cmp.FileNode):
                label = o.filename
            elif isinstance(o, cmp.FolderNode):
                label = o.foldername
            #elif isinstance(o, cmp.SharedInputNode):
            #    try:
            #        label = o.input_node.foldername
            #    except AttributeError:
            #        pass
            if label is None:
                label = "N/A"
            else:
                if o.input_node is not None:
                    label = format_path(label,
                                        path=path,
                                        replace_long_filename=True)
                else:
                    label = format_path(label)
            #if isinstance(o, cmp.FolderNode):
            #    o.gv_props["fillcolor"] = "red"
            s = node_str(o.id,
                         label=label,
                         **o.gv_props)
            return s, ""
        elif isinstance(o, cmp.SharedInputNode):
            label = "shared"
            o.gv_props["fillcolor"] = "slategray3"
            parent_scope = ""
            for node in o.output_nodes + [o.input_node]:
                if isinstance(node, (cmp.FileNode, cmp.FolderNode, cmp.SharedInputNode)):
                    s, p = traverse(node)
                    parent_scope += s
            s = node_str(o.id,
                         label=label,
                         **o.gv_props)
            return s, parent_scope
        elif isinstance(o, cmp.ComponentNode):
            label = o.get_name()
            if o.parallel:
                label += " (P)"
            if o.assoc_rna is not None:
                label += "\\n({})".format(o.assoc_rna)
            parent_scope = ""
            for node in o.output_nodes + [o.input_node]:
                if isinstance(node, (cmp.FileNode, cmp.FolderNode, cmp.SharedInputNode)):
                    s, p = traverse(node)
                    parent_scope += s
            s = node_str(o.id,
                         label=label,
                         **o.gv_props)
            return s, parent_scope
        else:
            print("Did not recognize class for {}".format(str(o)))
        return "", ""

    s, _ = traverse(pipeline)

    edges = pipeline.collect_edges()
    edges_str = ""
    for edge in edges:
        # shorten edges between non-intermediate FileNodes
        # and ComponentNodes
        intermediate = False
        if isinstance(edge.from_node, (cmp.FileNode, cmp.FolderNode)):
            if edge.from_node.input_node is not None:
                intermediate = True
        if isinstance(edge.to_node, (cmp.FileNode, cmp.FolderNode)):
            if len(edge.to_node.output_nodes)>0:
                intermediate = True
        if intermediate:
            #minlen = "2" # disabled for now
            minlen = "1"
        else:
            minlen = "1"

        # lengthen edges into certain components to help layout
        if isinstance(edge.to_node, cmp.ComponentNode):
            if ( any(edge.to_node.parent_component.get_name().startswith(n) for n in ["StarAligner",
                                                              "BowtieAligner",
                                                              "CalcProfile"]) and
                 "index" in edge.to_node.get_name() ):
                if not isinstance(edge.from_node, cmp.SharedInputNode):
                    minlen = "20"
            if isinstance(edge.to_node, cmp.SharedInputNode):
                minlen = "20"
            if edge.to_node.parent_component.get_name().startswith("CalcProfile"):
                minlen = "5"
            if edge.to_node.parent_component.get_name().startswith("RenderMutations"):
                minlen = "10"

        edges_str += edge_str(edge.from_node.id,
                              edge.to_node.id,
                              minlen=minlen)


    # add invisible edges between input and output nodes
    # to encourage natural top-to-bottom arrangement
    for c in pipeline.collect_low_level_components():
        for i in c.input_nodes:
            for j in c.output_nodes:
                edges_str += edge_str(i.id,
                                      j.id,
                                      style="invis",
                                      weight="100")

    s = body_str(legend_str(format_path(path)), s, edges_str)

    return s


def draw_flowchart(pipeline,
                   filename,
                   path="",
                   highlight_dir=None):

    p, ext = os.path.splitext(filename)

    dot = pipeline_to_dot(pipeline,
                          path=path,
                          highlight_dir=highlight_dir)

    # create path to output file if needed
    folder, name = os.path.split(filename)
    if len(folder) > 0:
        os.makedirs(folder, exist_ok=True)

    if False:
        # for debugging, output intermediate graphviz dot file
        dotpath = p + ".dot"
        o = open(dotpath, "w")
        o.write(dot)
        o.close()

    cmd = ["dot",
           "-T{}".format(ext[1:]),
           '>',
           "'" + filename + "'"]
    proc = sp.Popen(' '.join(cmd),
                    shell=True,
                    stdin=sp.PIPE,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    preexec_fn=os.setsid)
    stdout, stderr = proc.communicate(input=bytes(dot, 'UTF-8'))
    stdout = stdout.decode()
    stderr = stderr.decode()
    success = True
    if len(stdout.strip()) > 0:
        print("graphviz output:")
        print(stdout)
    err = stderr.strip()
    if len(err) > 0:
        # ignore font config warning (graphviz still performs render
        # even if FONTCONFIG_PATH not properly set).
        # also ignore weird library warning that doesn't seem to affect anything.
        if not err.startswith("Fontconfig error: Cannot load default config file") \
           and not err.endswith('libgvplugin_pango.so.6" - file not found'):
            print("graphviz error/warning message:")
            print(stderr.splitlines()[0])
            success = False
    return success
