#!/usr/bin/env python3
"""
Given a sequence name and a fasta file containing one or more sequences,
output the length of each sequence to separate files. Used to dynamically update
sequence length parameters for pipeline components such as MutationCounter.

"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import argparse, sys, os


def get_lengths(filename):
    f = open(filename, "r")
    lengths = []
    for line in f:
        if line[0] == ">":
            lengths.append(0)
        else:
            lengths[-1] += len(line.strip())
    return lengths


def write_lengths(outnames, lengths):
    assert len(outnames) == len(lengths)
    #o = os.path.split(prefix)[0]
    #if len(o)>0:
    #    os.makedirs(o, exist_ok=True)
    for i in range(len(lengths)):
        oname = outnames[i]
        o = os.path.split(oname)[0]
        if len(o)>0:
            os.makedirs(o, exist_ok=True)
        #oname = prefix+"l{}".format(i+1)
        o = open(oname, "w")
        o.write(str(lengths[i]))


if __name__=="__main__":
    parser = argparse.ArgumentParser()

    h = "Input fasta file"
    parser.add_argument("--fa", help=h, required=True, type=str)

    #h = 'Output prefix (will be followed by "l1", "l2", etc.)'
    #parser.add_argument("--prefix", help=h, required=False, type=str)

    h = 'Explicit output filenames (assumed same order as sequences in fasta)'
    parser.add_argument("--out", help=h, required=True, type=str, nargs="+")

    p = parser.parse_args(sys.argv[1:])

    lengths = get_lengths(p.fa)
    write_lengths(p.out, lengths)    
