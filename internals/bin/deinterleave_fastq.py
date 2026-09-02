#!/usr/bin/env python3
import sys, os
from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument("--input", type=str)
ap.add_argument('--unpaired-out', type=str)
ap.add_argument('--R1-out', type=str)
ap.add_argument('--R2-out', type=str)

# - assumes mixed paired/unpaired input
# - skips reads of length 1 (usually placeholders)


# FIXME: rewrite in c++
# FIXME: combine with deinterleave_fastq_columns.py as a more general utility

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


r1_out = open(pa.R1_out, 'w')
r2_out = open(pa.R2_out, 'w')
unpaired_out = open(pa.unpaired_out, 'w')


def gen():
    reads = []
    for r in iterate_fastq(pa.input):
        reads.append(r)
        if len(reads) == 2:
            r1_id = reads[0][0]
            r2_id = reads[1][0]
            if r2_id == r1_id:
                #reads[0][0] += ":R1"
                #reads[1][0] += ":R2"
                yield reads
                reads = []
            else:
                yield [reads[0]]
                reads = [reads[1]]
    if len(reads) > 0:
        yield reads

L = '@{}\n{}\n+\n{}\n'
for reads in gen():
    if len(reads) == 2:
        _, seq1, _ = reads[0]
        _, seq2, _ = reads[1]
        r1_out.write(L.format(*reads[0]))
        r2_out.write(L.format(*reads[1]))
        #if len(seq1) > 1 and len(seq2) > 1:
        #    r1_out.write(L.format(*reads[0]))
        #    r2_out.write(L.format(*reads[1]))
        #else:
        #    if len(seq1) > 1:
        #        unpaired_out.write(L.format(*reads[0]))
        #    if len(seq2) > 1:
        #        unpaired_out.write(L.format(*reads[1]))
    else:
        unpaired_out.writelines(L.format(*reads[0]))

r1_out.close()
r2_out.close()
unpaired_out.close()
