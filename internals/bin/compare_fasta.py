#!/usr/bin/env python3
"""
Check if two fasta files have identical sequences

"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys
import shutil
import traceback


def load_fasta(fastaname, rna=None, convert_to_rna=True):
    f = open(fastaname, "r")
    seq = ""
    rna_count = 0
    in_selected_rna = False
    seq_name = None
    if rna is None:
        in_selected_rna = True
    for line in f:
        if line.startswith(">"):
            name = line[1:].strip()
            rna_count += 1
            if len(seq) > 0:
                break
            if rna_count > 1 and rna is None:
                s = "Error: fasta file "
                s += "\"" + fastaname + "\""
                s += " contains more than one sequence, but no"
                s += " sequence name was specified."
                raise RuntimeError(s)
            if name == rna:
                in_selected_rna = True
            seq_name = name
            continue
        if in_selected_rna:
            seq += line.strip()
    if convert_to_rna:
        seq = seq.replace("T", "U")
    return seq_name, seq

try:
    if len(sys.argv)<3:
        print("Usage: python compare_fasta.py <seq1>.fa <seq2>.fa")
        sys.exit()

    _, s1 = load_fasta(sys.argv[1])
    _, s2 = load_fasta(sys.argv[2])

    term_width, term_height = shutil.get_terminal_size()
    s = "".join(['=' for x in range(term_width)])
    #print(s)
    if s1 == s2:
        #sys.stdout.write("PASS: Corrected sequence matches original sequence.\n")
        sys.exit(0)
    else:
        sys.stderr.write("ERROR: Corrected sequence does not match original sequence.\n")
        sys.stderr.write(s1+"\n")
        sys.stderr.write(s2+"\n")
        sys.exit(0)
    #print(s)
except Exception as e:
    if isinstance(e, KeyboardInterrupt):
        raise Exception(e)
    sys.stderr.write("ERROR: "+traceback.format_exc())
    sys.stderr.write("{}".format(e))
