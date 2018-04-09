"""
Render individual sequencing reads and mutations in their aligned
context, mostly for debugging purposes.
"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys
import operator


if len(sys.argv)<6:
    # FIXME: add some actual arg parsing
    print("Usage: python render_mutations_ps.py <aligned.sam> <target_name> <min_mapq> <parsed_mutations.mut> <classified_mutations.txt> <output.ps>")
    sys.exit()

# FIXME: add some labels to indicate what each row is showing

sam_file = open(sys.argv[1], "rU")
target_name = sys.argv[2]
min_mapq = int(sys.argv[3])
parsed_mutations_file = open(sys.argv[4], "rU")
classified_mutations_file = open(sys.argv[5], "rU")

o = open(sys.argv[6], "w")

char_height = 4
char_width = 3

vert_char_count = 20

starting_left = 75
starting_top = 750


class Mutation:
    def __init__(self, left, right, seq, qual, tag):
        self.left = int(left)
        self.right = int(right)
        self.seq = seq.replace('"','')
        self.qual = qual.replace('"','')
        self.tag = tag.replace('"','')


def mut_to_ps(line_number, title, left, right, read, parsed_mutations, classified_mutations):
    parsed_mutations = sorted(parsed_mutations, key=operator.attrgetter('left'))    
    classified_mutations = sorted(classified_mutations, key=operator.attrgetter('left'))  

    y_upper = starting_top - (line_number*vert_char_count*char_height)
    y_lower = y_upper - char_height
    y_lower2 = y_lower - char_height
    y_upper2 = y_upper + char_height
    y_upper_down = y_upper - char_height/2.0
    y_upper2_down = y_upper2 - char_height/2.0

    y_upper = int(round(y_upper))
    y_upper2 = int(round(y_upper2))
    y_upper2_down = int(round(y_upper2_down))
    y_upper_down = int(round(y_upper_down))
    y_lower = int(round(y_lower))
    y_lower2 = int(round(y_lower2))

    dim_colors = ["\n0 0 0.5 setrgbcolor",
                  "\n0 0.5 0 setrgbcolor",
                  "\n0.5 0 0 setrgbcolor",]
    bright_colors = ["\n0.2 0.2 1 setrgbcolor",
                     "\n0.2 1 0.2 setrgbcolor",
                     "\n1 0.2 0.2 setrgbcolor"]

    ref = "\n0 0 0 setrgbcolor"
    ref += "\n%i %i moveto\n(%s) show"%(starting_left-(len(title)+1)*char_width, int(round(y_lower-char_height/2.0)), title)
    mut = ""
    mut_bounds = ""

    for i in range(len(read)):
        x = int(round(starting_left + (i)*char_width))
        ref += "\n%i %i moveto\n(%s) show"%(int(round(x)), int(round(y_lower2)), read[i])

    c = 0
    for i in range(len(parsed_mutations)):
        mutation = parsed_mutations[i]
        if c >= len(bright_colors):
            c = 0
        mut += dim_colors[c]
        mut_bounds += bright_colors[c]        
        x = int(round(starting_left + (mutation.left+1-left)*char_width))
        x2 = int(round(starting_left + (mutation.right-left)*char_width))
        # rotate text if it overlaps adjacent mutation
        rot = 0
        try:
            next_left = parsed_mutations[i+1].left
            if mutation.left+len(mutation.seq) > next_left:
                rot = 45
        except IndexError:
            pass
        if x != x2:
            mut += "\n%i %i %i (%s) outputtext"%(int(round(x)), int(round(y_upper)), rot, mutation.seq)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x, y_lower,
                                                                  x, y_upper_down)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x, y_upper_down,
                                                                  x2, y_upper_down)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x2, y_lower,
                                                                  x2, y_upper_down)
        else:
            # draw slightly taller if a pure insert (otherwise can get hidden pretty easily)
            mut += "\n%i %i %i (%s) outputtext"%(int(round(x)), int(round(y_upper2)), rot, mutation.seq)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x, y_lower,
                                                                  x, y_upper2_down)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x, y_upper2_down,
                                                                  x2, y_upper2_down)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x2, y_lower,
                                                                  x2, y_upper2_down)
        c += 1

    # FIXME: change confusing var names to something generic
    y_upper = starting_top - (line_number*vert_char_count*char_height) - 3.5*char_height
    y4 = y_upper - char_height
    y5 = y_upper - 2*char_height
    y_lower = y_upper + char_height
    y_lower2 = y_lower + char_height
    y_upper2 = y_upper - char_height
    y6 = y_upper2 - char_height/2.0
    y7 = y6 - char_height/2.0
    y_upper_down = y_upper + char_height/2.0
    y_upper2_down = y_upper2 + char_height/2.0

    y_upper = int(round(y_upper))
    y4 = int(round(y4))
    y5 = int(round(y5))
    y6 = int(round(y6))
    y7 = int(round(y7))
    y_upper2 = int(round(y_upper2))
    y_upper2_down = int(round(y_upper2_down))
    y_upper_down = int(round(y_upper_down))
    y_lower = int(round(y_lower))
    y_lower2 = int(round(y_lower2))
    
    c = 0
    for i in range(len(classified_mutations)):
        mutation = classified_mutations[i]
        if c >= len(bright_colors):
            c = 0
        mut += dim_colors[c]
        mut_bounds += bright_colors[c]        
        x = int(round(starting_left + (mutation.left+1-left)*char_width))
        x2 = int(round(starting_left + (mutation.right-left)*char_width))
        # rotate text if it overlaps adjacent mutation
        rot = 0
        try:
            next_left = classified_mutations[i+1].left
            if mutation.left+len(mutation.seq) > next_left:
                rot = -45
        except IndexError:
            pass
        if x != x2:
            mut += "\n%i %i %i (%s) outputtext"%(x, y4, rot, mutation.seq)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x, y_lower,
                                                                  x, y_upper_down)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x, y_upper_down,
                                                                  x2, y_upper_down)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x2, y_lower,
                                                                  x2, y_upper_down)
            mut += "\n%i %i -45 (%s) outputtext"%(x, y5, mutation.tag)
        else:
            # draw slightly taller if a pure insert (otherwise can get hidden pretty easily)
            mut += "\n%i %i %i (%s) outputtext"%(x, y6, rot, mutation.seq)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x, y_lower,
                                                                  x, y_upper2_down)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x, y_upper2_down,
                                                                  x2, y_upper2_down)
            mut_bounds += "\n%i %i moveto\n%i %i lineto\nstroke"%(x2, y_lower,
                                                                  x2, y_upper2_down)
            mut += "\n%i %i -45 (%s) outputtext"%(x, y7, mutation.tag)        
        c += 1
    return [mut_bounds, ref, mut]

mut_classes = [
            "A-", "T-", "G-", "C-",
            "-A", "-T", "-G", "-C", "-N",
            "AT", "AG", "AC",
            "TA", "TG", "TC",
            "GA", "GT", "GC",
            "CA", "CT", "CG",
            "multinuc_deletion",
            "multinuc_insertion",
            "multinuc_mismatch",
            "complex_deletion",
            "complex_insertion",
            ]
new_classes = []
for c in mut_classes:
    new_classes += [c, c+"_ambig"]
mut_classes = new_classes
mut_classes += ["other"]
class_counts = {c:0 for c in mut_classes}
other_classes = []

header = """%!
%%DocumentMedia: custom 2000 842 80 () ()
/outputtext {
   /data exch def
   /rot exch def
   /y1 exch def
   /x1 exch def
   x1 y1 moveto
   rot rotate
   data show

   rot neg rotate
} def

1 1 1 setrgbcolor clippath fill
/mono findfont
5 scalefont setfont
"""

o.write(header+"\n")
o.write("\n%%%%Page: %i %i\n"%(1, 1))
start_line = 1
max_lines = 5000
last_page = 0
i = 0
total_nucs = 0
#for i in range(start_line+max_lines):
while True:
    cluster_id = "N/A"
    try:
        while True:
            l = sam_file.readline()
            if l[0] != '@':
                s = l.strip().split('\t')
                cluster_id = ":".join(s[0].split(":")[5:])
                cigarString = s[5]
                read = s[9]
                startIndex = int(s[3])-1
                target = s[2]
                mapq = int(s[4])
                if target == '*' or mapq < min_mapq or target != target_name:
                    continue
                break
    except IndexError:
        break
    #if cluster_id=="17761:3505":
    #    print("cluster {} is read {}".format(cluster_id, i+1))

    s = parsed_mutations_file.readline().strip().split()
    left = int(s[1])
    right = int(s[2])
    read = s[3]
    qual = s[4]
    parsed_mutations = []
    try:
        k = s[5:]
        while True:
            parsed_mutations.append(Mutation(k[0],k[1],k[2],k[3],''))
            del k[:4]
    except IndexError:
        pass

    s = classified_mutations_file.readline().strip().split()
    classified_mutations = []
    try:
        del s[:5] # first five fields are read name, effective read left and right
                  # target bounds, effective depth, and effective mutation count
        while True:
            classified_mutations.append(Mutation(s[0],s[1],s[2],s[3],s[4]))
            classif = s[4][1:-1]
            seq = s[2][1:-1]
            qual = s[3][1:-1]
            if classif in class_counts:
                class_counts[classif] += 1
            else:
                class_counts["other"] += 1
                if classif not in other_classes:
                    other_classes.append(classif)
            del s[:5]
    except IndexError:
        pass   

    ### debug
    if False:
        found_loc = False
        for m in classified_mutations:
            if m.right == 171:
                found_loc = True
                break
        if not found_loc:
            continue
    ###

    i += 1
    
    total_nucs += abs(right-left)+1

    lines_per_page = 700/(vert_char_count*char_height)

    if i < start_line or i > start_line+max_lines:
        continue

    page_num = int((i-start_line+1)/lines_per_page)
    line_number = i-start_line+1-page_num*lines_per_page

    if page_num > last_page:
        o.write("\nshowpage\n%%%%Page: %i %i\n"%(page_num+1, page_num+1))
    last_page = page_num

    title = "{}-{}".format(i, cluster_id)


    lines = mut_to_ps(line_number, title, left, right, read, parsed_mutations, classified_mutations)
    o.writelines(lines)


o.write("\nshowpage") # final showpage

print("Total reads: {}".format(i))
print("Total nts: {}".format(total_nucs))
print("Total classified mutation counts:")
for c in mut_classes:
    print("{}: {}".format(c, class_counts[c]))
print("Other mutation classes:")
for c in other_classes:
    print(str(c))

