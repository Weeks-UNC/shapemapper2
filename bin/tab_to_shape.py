"""
Given a tab-delimited reactivity summary file
output a .map file or .shape file
"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys
import os
import argparse

import numpy as np
from numpy import isfinite, nan


def load_tab(filename):
    f = open(filename, "rU")

    # do one pass to determine array length
    # TODO: might actually be faster to just resize array in memory and read in one pass
    length = 0
    f.readline()  # skip header
    for line in f:
        length += 1
    f.seek(0)

    if length == 0:
        s = "Error: file "
        s += "\"" + filename + "\""
        s += " contains no data."
        raise RuntimeError(s)

    line = f.readline()
    if line is None:
        raise RuntimeError("File \""+filename+"\" is empty.")
    headers = line.strip().split('\t')
    if len(headers)<1:
        raise RuntimeError("File \""+filename+"\" does not contain the expected header.")
    #expected_fields = ['HQ_profile','HQ_stderr', 'Norm_profile', 'Norm_stderr']
    # load normalized columns if present, otherwise load unnormed
    shape_header = None
    stderr_header = None
    seq_header = 'Sequence'
    if 'Norm_profile' in headers and 'Norm_stderr' in headers:
        shape_header = 'Norm_profile'
        stderr_header = 'Norm_stderr'
    elif 'HQ_profile' in headers and 'HQ_stderr' in headers:
        shape_header = 'HQ_profile'
        stderr_header = 'HQ_stderr'
    if shape_header is None or seq_header not in headers:
        raise RuntimeError("File \""+filename+"\" does not contain the expected columns.")
    profile_col = headers.index(shape_header)
    stderr_col = headers.index(stderr_header)
    seq_col = headers.index(seq_header)

    profile = np.empty(length)
    stderrs = np.empty(length)
    seq = np.empty(length, dtype=str)

    for i in range(length):
        line = f.readline()
        s = line.strip().split('\t')
        profile[i] = float(s[profile_col])
        stderrs[i] = float(s[stderr_col])
        seq[i] = s[seq_col].upper()

    return seq, profile, stderrs


def make_tree(path):
    # create path to file if needed
    d = os.path.split(path)[0]
    if len(d) > 0:
        os.makedirs(d, exist_ok=True)


def write_shape(shape, path):
    # write .SHAPE file
    make_tree(path)
    file_out = open(path, "w")
    nuc = 1
    for value in shape:
        if not isfinite(value):
            s = '-999'
        else:
            s = "{0:.6f}".format(value)
        file_out.write("{}\t{}\n".format(nuc, s))
        nuc += 1
    file_out.close()


def write_map(seq, shape, stderr, path):
    # write .map file
    make_tree(path)
    file_out = open(path, "w")
    nuc = 1
    for i in range(len(shape)):
        v1 = shape[i]
        v2 = stderr[i]
        if not isfinite(v1):
            s1 = -999
            s2 = 0
        else:
            s1 = "{0:.6f}".format(v1)
            s2 = "{0:.6f}".format(v2)
        file_out.write("{}\t{}\t{}\t{}\n".format(nuc, s1, s2, seq[i]))
        nuc += 1
    file_out.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    h = "Tab-delimited input file containing reactivity profile and stderrs."
    parser.add_argument("--infile", required=True, type=str, help=h)

    parser.add_argument("--map", required=False, type=str)
    parser.add_argument("--shape", required=False, type=str)

    p = parser.parse_args(sys.argv[1:])

    seq, shape, stderr = load_tab(p.infile)

    if p.map:
        write_map(seq, shape, stderr, p.map)
    if p.shape:
        write_shape(shape, p.shape)
