"""
Split a SAM file into multiple files by mapped target.

"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys, os
import argparse


# TODO: quantify number of paired reads that map to different RNAs (rough metric for PCR issues?)

def split_sam(sam=None,
              names=None,
              outs=None):
    if len(outs) != len(names):
        raise RuntimeError("Error: number of output files must match number of sequence target names.")
    f = open(sam, "rU")
    o = [open(x, "w") for x in outs]
    for line in f:
        # copy any headers/comment lines to all outputs
        if line[0] == '@':
            for x in o:
                x.write(line)
        else:
            try:
                target = line.split('\t')[2]
            except IndexError:
                raise RuntimeError("Error: SAM file appears misformatted")
            try:
                o[names.index(target)].write(line)
            except ValueError:
                # skip unmapped reads or reads mapped to sequences
                # not provided in names
                pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    h = "Input SAM file (uncompressed only)"
    parser.add_argument("-i", required=True, type=str, help=h)

    h = "List of sequence names"
    parser.add_argument("-n", required=True, type=str, nargs="*", help=h)

    h = "List of output files, in same order as sequence names"
    parser.add_argument("-o", required=True, type=str, nargs="*", help=h)

    p = parser.parse_args(sys.argv[1:])

    split_sam(sam=p.i,
              names=p.n,
              outs=p.o)