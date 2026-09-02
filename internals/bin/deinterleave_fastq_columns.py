#!/usr/bin/env python3
import sys, os

# FIXME: rewrite in c++

def iterate_fastq(filename):
    f = open(filename, "r")
    lines = []
    for line in f:
        # suppress jdb socket message
        if line.startswith("Listening"):
            continue
        lines.append(line)
        if len(lines) == 4:
            l = [lines[0][1:].split()[0], lines[1].strip(), lines[3].strip()]
            yield l
            lines = []

input = open(sys.argv[1], 'r')
r1_out = open(sys.argv[2], 'w')
r2_out = open(sys.argv[3], 'w')
unpaired_out = open(sys.argv[4], 'w')

for line in input:
    fields = line.rstrip().split('\t')
    # check for R1
    has_R1 = False
    try:
        if fields[0] != "" and len(fields[1]) > 1:
            has_R1 = True
    except IndexError:
        pass
    has_R2 = False
    try:
        if fields[4] != "" and len(fields[5]) > 1:
            has_R2 = True
    except IndexError:
        pass
    paired = (has_R1 and has_R2)

    if paired:
        r1_out.writelines(fields[0:4])
        r2_out.writelines(fields[4:8])
    else:
        if has_R1:
            unpaired_out.writelines(fields[0:4])
        elif has_R2:
            unpaired_out.writelines(fields[4:8])

