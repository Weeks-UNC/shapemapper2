#!/usr/bin/env python3
import sys, os
from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument("--input", type=str)
ap.add_argument("--R1", type=str)
ap.add_argument("--R2", type=str)
ap.add_argument('--output', type=str)

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

pa = ap.parse_args()

single_input = False
if pa.input:
    single_input = True

o = open(pa.output, "w", 1)

def gen():
    if single_input:
        reads = []
        for r in iterate_fastq(pa.input):
            if len(reads) == 0:
                reads.append(r)
            else:
                r1_id = reads[0][0]
                r2_id = r[0]
                if r2_id == r1_id:
                    reads.append(r)
                    #reads[0][0] += ":R1"
                    #reads[1][0] += ":R2"
                    yield reads
                    reads = []
                else:
                    yield reads
                    reads = [r]
        if len(reads) > 0:
            yield reads
    else:
        for r1, r2 in zip(iterate_fastq(pa.R1),
                          iterate_fastq(pa.R2)):
            #r1[0] += ":R1"
            #r2[0] += ":R2"
            yield [r1, r2]


for reads in gen():
    fields = []
    for r in reads:
        name, seq, qual = r
        fields += r
        #if len(seq) > 1:
        #    fields += r
    if len(fields) > 0:
        o.write('\t'.join(fields) + '\n')

o.close()
