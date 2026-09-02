#!/usr/bin/env python3
import sys, os
from itertools import zip_longest

# FIXME: rewrite in c++

def decode_flags(s):
    flags = [False]*12
    for i, c in enumerate(bin(int(s))[2:][::-1]):
        if c == '1':
            flags[i] = True
    return flags

# skip header lines, keep mate pairs together
def iterate_sam(f):
    lines = []
    for line in f:
        if line[0] != '@':
            flags = decode_flags(line.split()[1])
            paired = flags[0]
            first_in_pair = flags[6]
            second_in_pair = flags[7]
            mate_unmapped = flags[3]
            lines.append(line)
            if (not paired) or mate_unmapped:
                yield lines
                lines = []
            elif len(lines) == 2:
                yield lines
                lines = []
    if len(lines) > 0:
        yield lines

o = open(sys.argv[3], "w")
f1 = open(sys.argv[1], "r")
f2 = open(sys.argv[2], "r")

for lines1, lines2 in zip_longest(iterate_sam(f1),
                                  iterate_sam(f2)):
    if lines1 is not None:
        o.writelines(lines1)
    if lines2 is not None:
        o.writelines(lines2)

