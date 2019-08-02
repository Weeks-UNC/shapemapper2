#!/usr/bin/env python3
import sys, os

# FIXME: rewrite in c++

def iterate_fastq(filename):
    f = open(filename, "rU")
    lines = []
    for line in f:
        lines.append(line.strip())
        if len(lines) == 4:
            yield lines
            lines = []

o = open(sys.argv[3], "w")

for r1, r2 in zip(iterate_fastq(sys.argv[1]),
                  iterate_fastq(sys.argv[2])):
    #r1[0] += ":R1"
    #r2[0] += ":R2"
    o.write('\n'.join(r1+r2) + '\n')

