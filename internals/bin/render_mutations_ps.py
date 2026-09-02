#!/usr/bin/env python3
"""
Render individual sequencing reads and mutations in their aligned
context, mostly for debugging purposes.
"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
import sys
import string
import operator
from argparse import ArgumentParser as AP

# "[read]"
# optional tag or string indicating type of read (UNPAIRED, UNPAIRED_R1, UNPAIRED_R2, PAIRED_R1, PAIRED_R2, MERGED)
# left mapping position (0-based)
# right mapping position (0-based)
# optional strand
# reference sequence
# optional quality scores aligned to reference sequence
# optional depth (local_effective_depth)
# optional inferred adduct locations (local_effective_count)
# optional associated mutations and quality scores

# "[separator] " - visual separator

# any non-bracketed string - string indicating computation just performed or other info
#       like cluster ID

# concordant R1 and R2 are processed in pairs, and rendered
# so right-most read is offset to line up with left-most read reference

# FIXME: add a more complete legend at the beginning

show_bounding_box = False
use_full_qual_gradient = False
only_pairs = False

# FIXME: make a param class
def init_globals(maxins=800, mustspan=None):
    global page_index
    global primers
    global char_height
    global char_width
    global row_width
    global margin
    global page_width
    global page_height
    global max_height
    global starting_top
    global starting_left
    global max_pages
    global vert_char_count
    global row_height
    global rows_per_page
    global must_span

    page_index = 0
    primers = []

    # all postscript units in points: 1/72 inch
    char_height = 4
    char_width = 3

    #a4_aspect = 595. / 842.
    a4_aspect = 842. / 595. # landscape orientation
    max_insert = maxins
    row_width = max_insert * char_width
    margin = 75
    page_width = row_width + 2*margin
    page_height = float(page_width) / a4_aspect
    max_height = page_height - 2*margin
    starting_top = page_height - margin
    starting_left = margin
    max_pages = 100
    vert_char_count = 9
    row_height = vert_char_count*char_height
    # adjust rows_per_page to maintain roughly A4 aspect ratio
    rows_per_page = max_height/(vert_char_count*char_height)

    if mustspan is not None:
        mustspan = list(map(int, mustspan.replace(',',' ').replace('-',' ').split()))
        assert len(mustspan) == 2
    must_span = mustspan


block_count = 0

class Mutation:
    def __init__(self, left, right, seq, qual, tag):
        self.left = int(left)
        self.right = int(right)
        self.seq = seq.replace('"','')
        self.qual = qual.replace('"','')
        self.tag = tag.replace('"','')


def parse_mutation_fields(field):
    mutations = []
    if len(field) == 0:
        return mutations
    s = field.split(' ')
    while len(s) > 0:
        mutations.append(Mutation(s[0], s[1], s[2], s[3], s[4]))
        del s[:5]
    return mutations


class Read:
    def __init__(self, line):
        s = line.rstrip('\n').split('\t')
        self.read_type, \
        self.left, \
        self.right, \
        self.strand, \
        self.mapping_category, \
        self.primer_pair, \
        self.seq, \
        self.qual, \
        self.mapped_depth, \
        self.depth, \
        self.count = s[1:12]
        self.left = int(self.left)
        self.right = int(self.right)
        self.primer_pair = int(self.primer_pair)
        self.mutations = parse_mutation_fields(s[12])

        # quick hack so mapped depths get shown for initial
        # processing stages
        if self.depth == "":
            self.depth = self.mapped_depth

        #print(self.read_type)


class Primer:
    def __init__(self):
        self.seq = ""
        self.strand = '+'
        self.left = 0
        self.right = 0

    def __str__(self):
        return "{} {} {} {}".format(self.seq, self.strand, self.left, self.right)

class PrimerPair:
    def __init__(self):
        self.fw = Primer()
        self.fw.strand = '+'
        self.rv = Primer()
        self.rv.strand = '-'

    def __str__(self):
        return "pair:\n{}\n{}\n".format(self.fw, self.rv)


def load_primers(filename):
    if filename is None or filename == "":
        return []

    f = open(filename, "r")
    primers = []
    for line in f:
        if line[0] == '>' or len(line) < 1:
            continue
        pair = PrimerPair()
        s = line.rstrip().split()
        if line[0] in 'AUTGCN':
            pair.fw.seq = s[0]
            pair.rv.seq = s[1][::-1] # store reverse primer sequence reversed to simplify rendering
            primers.append(pair)
        else:
            primers[-1].fw.left = int(s[0])
            primers[-1].fw.right = int(s[1])
            primers[-1].rv.left = int(s[2])
            primers[-1].rv.right = int(s[3])
    return primers



def set_color(*args):
    if len(args) == 3:
        r, g, b = args
    else:
        r, g, b = args[0]
    color_template = "\n{:.3f} {:.3f} {:.3f} setrgbcolor"
    return color_template.format(r,g,b)


def save_color():
    return "\n/savergbcolor {{ currentrgbcolor }} def"


def restore_color():
    return "\nsavergbcolor setrgbcolor"

def save_state():
    return "\ngsave"

def restore_state():
    return "\ngrestore"

# FIXME: maybe use scale command to allow object placement with > pt precision
#        - smaller graphics and the ability to fit everything on standard size page
# - for now just using pt sizes but at least downscaling afterward to fit to pdf page
def translate(loc):
    left, top = loc
    return "\n{left:.0f} {top:.0f} translate".format(**locals())

def text(x=0, y=0, rot=0, s='', color=None):
    text_template = ""
    if color is not None:
        text_template += save_color()
        text_template += set_color(color)
    text_template += "\n{x:.0f} {y:.0f} {rot:.0f} ({s}) outputtext"
    if color is not None:
        text_template += restore_color()
    return text_template.format(**locals())


def line_stroke(x1=0, y1=0, x2=0, y2=0, lw=1, color=None):
    line_template = ""
    if lw != 1:
        line_template += "\n/savelinewidth currentlinewidth def"
        line_template += "\n{lw:.0f} setlinewidth"
    if color is not None:
        line_template += save_color()
        line_template += set_color(color[0], color[1], color[2])
    line_template += "\n{x1:.0f} {y1:.0f} moveto\n{x2:.0f} {y2:.0f} lineto\nstroke"
    if lw != 1:
        line_template += "\nsavelinewidth setlinewidth"
    if color is not None:
        line_template += restore_color()
    return line_template.format(**locals())


def rect(x=0, y=0, w=0, h=0, filled=True, color=None, **kwargs):
    template = ""
    if color is not None:
        template += save_color()
        template += set_color(color)
    template += "\nnewpath {x:.0f} {y:.0f} {w:.0f} {h:.0f} "
    if filled:
        template += "rectfill"
    else:
        template += "rectstroke"
    if color is not None:
        template += restore_color()
    return template.format(**locals())


def triup(x=0, y=0, w=0, h=0, **kwargs):
    template = "\nnewpath {x:.0f} {y:.0f} {w:.0f} {h:.0f} triup"
    return template.format(**locals())

def triright(x=0, y=0, w=0, h=0, **kwargs):
    template = "\nnewpath {x:.0f} {y:.0f} {w:.0f} {h:.0f} triright"
    return template.format(**locals())

def trileft(x=0, y=0, w=0, h=0, **kwargs):
    template = "\nnewpath {x:.0f} {y:.0f} {w:.0f} {h:.0f} trileft"
    return template.format(**locals())


def draw_separator(color=None):
    y = char_height/2.
    height = 2*char_height
    s = "\n%SEPARATOR%\n"
    s += line_stroke(x1=-10, y1=y, x2=row_width, y2=y, lw=2, color=color)
    if show_bounding_box:
        s += text(x=0, y=0, s="*", color=red)
        s += rect(x=0, y=0, w=row_width, h=height, filled=False, color=magenta)

    return s, height

def draw_na():
    return draw_separator(color=(0.8,0.5,0.8))

red = (1,0,0)
magenta = (1,0,1)


def draw_comment(s):
    rendered = text(x=0, y=0, s=s)
    h = 2*char_height
    if show_bounding_box:
        rendered += text(x=0, y=0, s="*", color=red)
        rendered += rect(x=0, y=0, w=row_width, h=h, filled=False, color=magenta)

    return rendered, h


lo_qual = 15.0 # 15 pretty reasonable
hi_qual = 40.0
gradient_name = "magma"
use_seaborn_palette = False

if use_full_qual_gradient:
    from matplotlib import pyplot as plt
    import matplotlib.colors as mp_colors
    import matplotlib.cm as mp_cm
    mp_gradients = [
        'viridis', 'plasma', 'inferno', 'magma',

        'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',

        'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
        'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
        'hot', 'afmhot', 'gist_heat', 'copper',
    ]
    if use_seaborn_palette:
        # note: seaborn not currently installed in thirdparty conda packages
        from seaborn import cubehelix_palette as cubehelix
        cmap = cubehelix_palette(n_colors=int(round(hi_qual-lo_qual+1)),
                                 start=0,
                                 rot=1.0,
                                 hue=0.8,
                                 light=0.9,
                                 dark=0.1,
                                 reverse = True,
                                 as_cmap=True)
    else:
        cmap = plt.get_cmap(gradient_name)
    cnorm = mp_colors.Normalize(vmin=lo_qual, vmax=hi_qual)
    scalar_map = mp_cm.ScalarMappable(norm=cnorm, cmap=cmap)
    def grad_color(val):
        if (val == 0.0) or (val > 40.0):
            return set_color(1,1,1)
        r,g,b,a = scalar_map.to_rgba(val)
        return set_color(r,g,b)
else:
    # setup a gradient from blue-white at qual 40 to
    # some sort of gray-red at qual 0 (excluding exact qual 0,
    # used as a placeholder in reference sequence if no basecall present, as in a deletion)
    lo_color = (0.6, 0.2, 0.2)
    hi_color = (0.9, 0.95, 1.0)

    def grad_color(val):
        rgb = []
        if (val == 0.0) or (val > 40.0):
            return set_color(1,1,1)
        adj_val = min(max(val, lo_qual), hi_qual)
        frac = (adj_val-lo_qual)/(hi_qual-lo_qual)
        for i in range(3):
            comp = lo_color[i] + frac*(hi_color[i]-lo_color[i])
            rgb.append(comp)
        return set_color(*rgb)

def draw_quals(quals, rot=0):
    s = ""

    # - for text(), rotation is handled in the custom PS routine "outputtext"
    # - probably easier to handle color gradient here rather than in PS routine
    if rot != 0:
        s += "\n/rot {rot} def".format(rot=rot)
    for i in range(len(quals)):
        phred = ord(quals[i])-33
        x = i * char_width
        w = char_width
        y = 0
        h = char_height
        if rot != 0:
            s += "\nrot rotate"
        s += grad_color(phred)
        s += rect(**locals())
        if rot != 0:
            s += "\nrot neg rotate"
    s += set_color(0,0,0)
    return s




def draw_sequence(seq, qual=None):
    s = ""

    # FIXME: offset is smaller than a point
    if qual is not None:
        #s += translate((0, -0.1*char_height))
        s += draw_quals(qual)
        #s += translate((0, 0.1*char_height))

    for i in range(len(seq)):
        x = i * char_width
        s += text(x=x,
                  y=0,
                  s=seq[i])
    return s


def draw_mutations(mutations, left, color_index=None):
    mut_qual_s = ""
    mut_s = ""
    mut_bounds_s = ""
    mut_type_s = ""
    bright_colors = [(0, 0, 0.5),
                     (0.2, 0.5, 0.3),
                     (0.5, 0, 0)]
    dim_colors = [(0,0,0.25),
                  (0.1, 0.25, 0.15),
                  (0.25, 0, 0)]

    # FIXME: any easy way to maintain more consistent colors between mate pairs?
    if color_index is None:
        color_index = 0
        try:
            color_index = mutations[0].left % 2
        except IndexError:
            pass

    for i in range(len(mutations)):
        mutation = mutations[i]

        mut_s += set_color(dim_colors[color_index])
        mut_bounds_s += set_color(bright_colors[color_index])

        x = (mutation.left+1-left)*char_width
        x2 = (mutation.right-left)*char_width
        # rotate text if it overlaps adjacent mutation
        rot = 0
        try:
            next_left = mutations[i+1].left
            if mutation.left+len(mutation.seq) > next_left:
                rot = 45
        except IndexError:
            pass
        seq = mutation.seq
        #if seq == "":
        #    seq = "-"
        ytext = 0
        yhigh = char_height
        if x == x2:
            # draw slightly taller if a pure insert (otherwise can get hidden pretty easily)
            ytext = 1.4*char_height
            yhigh = 2.4*char_height

        mut_qual_s += translate((x, ytext))
        mut_qual_s += draw_quals(mutation.qual, rot=rot)
        mut_qual_s += translate((-x, -ytext))

        mut_s += text(x=x,
                      y=ytext,
                      rot=rot,
                      s=seq)

        # FIXME: place more intelligently, maybe abbreviate and add key
        if mutation.tag != "":
            mut_type_s += text(x=x, y=yhigh, rot=20, s=mutation.tag, color=(0.5,0.5,0.5))

        # FIXME: this is smaller than a point, so gets rounded away. need to implement scaling to get more precision
        offset = char_width*0.01
        mut_bounds_s += line_stroke(x1=x+offset, y1=0,
                                    x2=x+offset, y2=yhigh)
        mut_bounds_s += line_stroke(x1=x+offset, y1=yhigh,
                                    x2=x2-offset, y2=yhigh)
        mut_bounds_s += line_stroke(x1=x2-offset, y1=0,
                                    x2=x2-offset, y2=yhigh)
        color_index += 1
        color_index = color_index if color_index < 3 else 0
    s = ""
    s += translate((0, char_height))
    s += mut_bounds_s
    s += translate((0, 1.2*char_height))
    s += mut_qual_s
    s += mut_s
    s += mut_type_s
    s += translate((0, -1.2*char_height))
    s += translate((0, -char_height))
    s += set_color(0,0,0)
    return s


def draw_depth(depth):
    s = ""
    for i in range(len(depth)):
        if depth[i] == '0':
            continue
        x = i * char_width
        w = char_width
        y = 0
        h = char_height
        s += rect(**locals())
    return s


def draw_counts(counts):
    # render inferred adduct locations as single triangles
    s = ""
    s += set_color(0.2,0.5,1.0)
    for i in range(len(counts)):
        if counts[i] == '0':
            continue
        x = (i + 0.5) * char_width
        w = char_width * 1.8
        y = 0
        h = char_height*1.8
        s += triup(**locals())
    s += set_color(0,0,0)
    return s


def draw_primers(read):
    s = ""
    # WIP: - add arrow indicating strand
    #        - handle overlaps, stack pairs if more than one within view
    #        - indicate primer pairs, maybe draw a light gray line connecting the primers
    inview_pairs = []
    if read.primer_pair >= 0:
        for i, p in enumerate(primers):
            if read.primer_pair == i:
                inview_pairs.append(p)
    else:
        for p in primers:
            if (((read.left <= p.fw.left <= read.right) or
                (read.left <= p.fw.right <= read.right)) or
               ((read.left <= p.rv.left <= read.right) or
                (read.left <= p.rv.right <= read.right))):
                inview_pairs.append(p)

    if len(inview_pairs) == 0:
        return "", 0.

    height = 0.8
    o = 1.0
    s += translate((0, 1.2*char_height))
    for n, pair in enumerate(inview_pairs):
        # line connecting forward and reverse primers in pair
        x1 = max((pair.fw.right + 2.5 - read.left), 0) * char_width
        x2 = min((pair.rv.left - 1.5 - read.left), read.right - read.left + 1) * char_width
        y1 = o * (n + 0.5) * char_height
        y2 = y1
        lw = 1.0
        color = (0.8, 0.8, 0.8)
        s += line_stroke(x1=x1, y1=y1, x2=x2, y2=y2, lw=lw, color=color)

        for k, p in enumerate([pair.fw, pair.rv]):
            color = (0.4,0,0)
            yo = 0.3
            if p.strand == '-':
                color = (0,0,0.4)
                yo = -0.3
            s += set_color(color)
            if ((read.left <= p.left <= read.right) or
                (read.left <= p.right <= read.right)):
                for i, c in enumerate(p.seq):
                    x = (p.left - read.left + i) * char_width
                    s += text(x=x,
                              y=o*(n+yo)*char_height,
                              s=c)
                y = o * (n + yo + 1) * char_height
                # k==0 is fw, 1 is rv
                if k == 0:
                    x = (p.right - read.left + 1.4) * char_width
                    s += triright(x=x, y=y, w=char_width, h=char_height)
                elif k == 1:
                    x = (p.left - read.left - 0.4) * char_width
                    s += trileft(x=x, y=y, w=char_width, h=char_height)
            s += set_color(0, 0, 0)

        s += set_color(0,0,0)

        height += o
    s += translate((0, -1.2*char_height))
    return s, height


def draw_read(read, color_index=None):
    s = ""
    # WIP: translate so bottom is at ~0

    s += draw_sequence(read.seq, qual=read.qual)
    s += draw_mutations(read.mutations, read.left)

    # FIXME: only draw primers once for read pairs
    sp, hp = draw_primers(read)
    if len(primers) > 0:
        s += translate((0, -(hp+1.2) * char_height))
        s += sp
        s += translate((0, (hp+1.2) * char_height))

    s += translate((0, -(1.2+hp)*char_height))
    s += draw_depth(read.depth)
    s += translate((0, (1.2+hp)*char_height))

    s += translate((0, -(1.4+hp)*char_height))
    s += draw_counts(read.count)
    s += translate((0, (1.4+hp)*char_height))

    t = "mapped {}-{}".format(read.left+1, read.right+1)
    if read.strand != "N/A":
        t += " ({})".format(read.strand)
    if read.read_type != "" and read.read_type != "UNSPECIFIED_READ_TYPE":
        t = read.read_type + ", " + t
    s += translate((0, 3*char_height))
    s += text(x=0,y=0, s=t)
    s += translate((0, -3*char_height))

    yo = 3.5+hp
    s = translate((4 * char_width, yo * char_height)) \
        + s \
        + translate((-4 * char_width, -yo*char_height))

    height = (8+hp)*char_height

    if show_bounding_box:
        s += text(x=0, y=0, s="*", color=red)
        s += rect(x=0, y=0, w=row_width, h=height, filled=False, color=magenta)

    return s, height






def draw_read_pair(reads):
    s = ""
    s1, h1 = draw_read(reads[0])
    s2, h2 = draw_read(reads[1])

    left = min(reads[0].left, reads[1].left)

    offset1 = (reads[0].left-left)*char_width
    offset2 = (reads[1].left-left)*char_width
    s += translate((offset1,h2))
    s += s1
    s += translate((-offset1,-h2))
    s += translate((offset2, 0))
    s += s2
    s += translate((-offset2, 0))

    height = h1+h2

    if show_bounding_box:
        s += text(x=0, y=0, s="*", color=red)
        s += rect(x=0, y=0, w=row_width, h=height, filled=False, color=magenta)

    return s, height



class Row:
    def __init__(self,
                 data,
                 row_type):
        self.data = data
        self.row_type = row_type
        assert row_type in ["Read", "ReadPair", "Comment", "Separator"]

    def render(self):
        if self.row_type == "Read":
            return draw_read(self.data)
        elif self.row_type == "ReadPair":
            return draw_read_pair(self.data)
        elif self.row_type == "Comment":
            return draw_comment(self.data)
        elif self.row_type == "Separator":
            return draw_separator()


def deindent(s, level, tab_size=4):
    lines = s.splitlines()
    wide_space = ''.join([' ']*tab_size)
    for i in range(len(lines)):
        for n in range(level):
            lines[i] = lines[i].lstrip(wide_space)
    return '\n'.join(lines)


def header():
    header = """%!
                %%DocumentMedia: custom {page_width} {page_height} 80 white ()
                /outputtext {{
                   /data exch def
                   /rot exch def
                   /y1 exch def
                   /x1 exch def
                   x1 y1 moveto
                   rot rotate
                   data show
                   rot neg rotate

                }} def

                /triup {{
                   /h exch def
                   /w exch def
                   /y exch def
                   /x exch def
                   % x,y should be the loc of the top pt of the triangle
                   x y moveto
                   /halfw w 2.0 div def
                   /x2 x halfw add def
                   /y2 y h sub def
                   x2 y2 lineto
                   /x3 x halfw sub def
                   /y3 y2 def
                   x3 y3 lineto
                   closepath
                   fill
                }} def

                /triright {{
                   /h exch def
                   /w exch def
                   /y exch def
                   /x exch def
                   % x,y should be the loc of the top left pt of the triangle
                   x y moveto
                   /halfh h 2.0 div def
                   /x2 x w add def
                   /y2 y halfh sub def
                   x2 y2 lineto
                   /x3 x def
                   /y3 y h sub def
                   x3 y3 lineto
                   closepath
                   fill
                }} def

                /trileft {{
                   /h exch def
                   /w exch def
                   /y exch def
                   /x exch def
                   % x,y should be the loc of the top right pt of the triangle
                   x y moveto
                   /halfh h 2.0 div def
                   /x2 x w sub def
                   /y2 y halfh sub def
                   x2 y2 lineto
                   /x3 x def
                   /y3 y h sub def
                   x3 y3 lineto
                   closepath
                   fill
                }} def

                0 0 0 setrgbcolor
                /Helvetica-Bold findfont
                5 scalefont
                setfont
            """
    header = header.format(page_height=page_height,
                           page_width=page_width)
    header = deindent(header, 4)

    return header

# FIXME: for page size stuff, see discussion at https://forums.adobe.com/thread/378001
# this doesn't seem to have any real effect, at least on ps2pdf conversion
def start_page(page_index):
    s = ""
    s += "\n\n%%Page: {page} {page}\n"
    s += "\n%%PageMedia: custom {page_width} {page_height} 80 white ()"
    s += "\n%%PageBoundingBox: {page_width} {page_height} dx dy\n"
    s = s.format(page=page_index+1,
                 page_width=page_width,
                 page_height=page_height)
    return s

def end_page():
    return "\nshowpage"


def parse_debug_file(f):
    global block_count
    n = 1
    for line in f:
        #print("line {}".format(n))
        #print(line.rstrip())
        #sys.stdout.flush()
        #if n % 100 == 0:
        #    print("line {}".format(n))
        #    print("block_count {}".format(block_count))
        if line.startswith('[read]'):
            yield Read(line)
        else:
            yield line.rstrip()
        n += 1


def group_objs(f):
    group = []
    for obj in parse_debug_file(f):
        if (isinstance(obj, Read) and
            (obj.read_type == "PAIRED_R1" or
             obj.read_type == "PAIRED_R2")):
            if len(group) > 0:
                if isinstance(group[0], Read):
                    group.append(obj)
                else:
                    yield group
                    group = [obj]
        else:
            if len(group) > 0:
                yield group
                group = []
            group = [obj]
        if len(group) == 2:
            yield group
            group = []
    if len(group) > 0:
        yield group


def classify_objs(f):
    block = []
    for group in group_objs(f):
        item = group[0]
        if len(group) == 2:
            block.append(Row(group, "ReadPair"))
        elif isinstance(item, Read):
            block.append(Row(item, "Read"))
        elif isinstance(item, str):
            if item.startswith('[separator]'):
                yield block
                block = [Row("", "Separator")]
            else:
                block.append(Row(item, "Comment"))
    if len(block) > 0:
        yield block


def filter_pairs(f, only_pairs=False):
    for block in classify_objs(f):
        #print('-#-'.join([row.row_type for row in block]))
        s = []
        for row in block:
            s.append("({})-({})")
            if row.row_type == "Read":
                s[-1] = s[-1].format(row.row_type, row.data.read_type)
            elif row.row_type == "ReadPair":
                s[-1] = s[-1].format(row.row_type, "{}:{}".format(row.data[0].read_type, row.data[1].read_type))
            else:
                s[-1] = s[-1].format(row.row_type, str(row.data))
        if only_pairs:
            found_pair = False
            for row in block:
                if row.row_type == "ReadPair":
                    found_pair = True
                    break
            if not found_pair:
                continue
            yield block
        else:
            yield block


def filter_span(f):
    global must_span

    for block in f:
        if must_span is None:
            yield block
        else:
            span = False
            for row in block:
                if row.row_type == "Read":
                    if (row.data.left <= must_span[0]) and (row.data.right >= must_span[1]):
                        span = True
                elif row.row_type == "ReadPair":
                    if ((min(row.data[0].left,row.data[1].left) <= must_span[0]) and
                        (max(row.data[0].right,row.data[1].right) >= must_span[1])):
                        span = True
                if span:
                    break
            if span:
                yield block
            else:
                continue


def render_rows(f):
    global block_count
    block_count = 0
    for block in filter_span(filter_pairs(f, only_pairs=only_pairs)):
        block_count += 1
        for row in block:
            yield row.render()


def legend():
    global page_index
    s = ""

    s += save_state()
    s += translate((starting_left, page_height*0.85))
    #quals = [chr(q+33) for q in list(range(0,41))]
    quals = list(range(0,41))

    s += '/Helvetica-Bold findfont\n7 scalefont\nsetfont\n'
    s += text(x=15, y=-20, rot=0, s='Basecall quality (phred score)')
    s += '/Helvetica-Bold findfont\n5 scalefont\nsetfont\n'

    for i in range(len(quals)):
        #phred = ord(quals[i])-33
        phred = quals[i]
        x = i * char_width
        w = char_width*1.2
        y = 0
        h = char_height
        s += grad_color(phred)
        s += rect(**locals())
        # ticks
        if phred % 5 == 0:
            s += set_color(0,0,0)
            s += line_stroke(x1=x+char_width/2.0, y1=y, x2=x+char_width/2.0, y2=y-char_height, lw=1)
            t = str(phred)
            rot=0
            y = y-char_height*2-1
            x = x-char_width*0.3
            if phred == 0:
                t = "none"
                rot = -90
                y += char_height
                x += char_width/2.0
            if phred == 0:
                s += '/Helvetica-Bold findfont\n4 scalefont\nsetfont\n'
            s += text(x=x, y=y, rot=rot, s=t)
            if phred == 0:
                s += '/Helvetica-Bold findfont\n5 scalefont\nsetfont\n'
    s += set_color(0,0,0)

    s += restore_state()
    s += end_page()
    page_index += 1
    s += start_page(page_index)

    return s



def main():
    global primers
    global max_pages
    global page_index

    ap = AP()
    ap.add_argument("--max-length", type=int, default=500)
    ap.add_argument("--input", type=str, required=True)
    ap.add_argument("--output", type=str, required=True)
    ap.add_argument("--primers", type=str, required=False, default="")
    ap.add_argument("--max-pages", type=int, default=100)
    ap.add_argument("--must-span", type=str, default=None)
    # FIXME: span argument probably isn't convenient for runs with multiple RNA targets
    p = ap.parse_args()

    init_globals(maxins=p.max_length,
                 mustspan=p.must_span)
    max_pages = p.max_pages
    primers = load_primers(p.primers)
    input_file = open(p.input, 'r')
    o = open(p.output, "w")

    o.write(header())
    start_loc = [starting_left, starting_top]
    cursor_loc = start_loc[:]
    page_index = 0
    o.write(start_page(page_index))
    o.write(legend())
    page_index += 1
    for s, h in render_rows(input_file):
        if starting_top-cursor_loc[1]-h > max_height:
            o.write(end_page())
            page_index += 1
            o.write(start_page(page_index))
            if page_index + 1 > max_pages:
                # consume the rest of the input stream so no early
                # terminations due to broken pipe
                #for line in input_file:
                #    pass
                # read 128kB at a time until EOF (read() returns empty string)
                while len(input_file.read(2**17)) > 0:
                    pass
                break
            cursor_loc = start_loc[:]
        o.write(save_state())
        o.write(translate(cursor_loc))
        # move row so upper left corner is at 0,0
        # (vs. the lower left corner at 0,0 as it within each draw function)
        s = translate((0,-h)) + s + translate((0,h))
        o.write(s)
        o.write(restore_state())
        o.flush()
        cursor_loc[1] -= h
    o.write(end_page()) # this will sometimes result in an extra blank page, but who cares
    o.write('\n%')
    o.flush()
    o.close()
    input_file.close()

if False:
    open("tmp.txt", "w").write(""">RNA-A
CTGGGACTTCCGAGGCAAC CATCACCTAGGAGGACGTACA
14 32 209 229
TGGGAAGGAGAGCGTCGTTA CAGTTCCAGGTGTCCTGCTT
147 166 336 355
GTCTGGTGGTGGGTCGTAAG GACAGTCGCTCCGTGACAG
419 438 593 611""")
    primers = load_primers("tmp.txt")
    for p in primers:
        print(p)
else:
    # run normally
    main()
    #os.system('gnome-open {1} && gedit {1}'.format(sys.argv[2]))
