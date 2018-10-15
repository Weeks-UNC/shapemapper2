# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys, os
from copy import deepcopy
from argparse import ArgumentParser
from util import sanitize

class Primer:
    def __init__(self,
                 seq,
                 left,
                 right):
        self.seq = seq
        self.left = left
        self.right = right


def complement(seq):
    d = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }
    comp = ""
    for c in seq:
        if c in d:
            comp += d[c]
        else:
            comp += 'N'
    return comp


def reverse_complement(seq):
    return complement(seq)[::-1]


def locate_primer_pairs(fastas,
                        rna_names,
                        primer_filenames,
                        primers_in_sequence):
    # Use primer sequences to find expected mapping locations for paired
    # amplicon reads
    # - If no RNA name provided, apply same primers to all target RNAs
    # - If primers_in_sequence, use lowercase sequence on either end of target
    #   sequence(s) as primer locations
    # - Error if RNA name not found in primers
    # - Error if RNA name in primers not found in targets
    # - Error if exact primer match with target sequence not found
    # - Error if multiple copies of primer pairs present in input
    assert isinstance(fastas, list)
    primers = {name:[] for name in rna_names+["__ALL_TARGETS__"]}  # dict of lists of pairs by RNA name

    current_RNA = "__ALL_TARGETS__"
    if primer_filenames is not None:
        for filename in primer_filenames:
            for line in open(filename, "rU"):
                if line[0] == '>':
                    current_RNA = sanitize(line[1:].rstrip())
                else:
                    s = line.strip().split()
                    if len(s) == 0:
                        continue
                    if len(s) != 2:
                        raise RuntimeError("Error: primer sequences must be specified with one pair per line, separated by spaces. The RNA name for a group of primers should be on the preceding line beginning with '>'.")
                    pair = (Primer(s[0],None, None), Primer(s[1], None, None))
                    if current_RNA not in primers:
                        raise RuntimeError("Error: RNA name '{}' in primer file does not match any RNA name provided in fasta file(s).".format(current_RNA))
                    primers[current_RNA].append(pair)

    # raise an error if duplicated primer sequences present
    for rna in primers:
        cat_seqs = []
        for p in primers[rna]:
            fw = p[0]
            rv = p[1]
            cat_seqs.append(fw.seq+rv.seq)
        if len(cat_seqs) > len(list(set(cat_seqs))):
            raise RuntimeError("Error: duplicated amplicon primer pair(s) present.")

    def iterate_fasta(filenames):
        for filename in filenames:
            f = open(filename, "rU")
            name = None
            seq = ""
            for line in f:
                if line[0] == '>':
                    if len(seq) > 0:
                        yield name, seq
                    name = sanitize(line[1:].rstrip())
                    seq = ''
                else:
                    seq += ''.join(line.strip().split())
            if len(seq) > 0:
                yield name, seq

    # identify primers from lowercase sequence on either end of each RNA
    if primers_in_sequence:
        for name, seq in iterate_fasta(fastas):
            fw_seq = ''
            i = 0
            while i < len(seq):
                c = seq[i]
                if c.islower():
                    fw_seq += c
                else:
                    break
                i += 1
            fw_primer = Primer(fw_seq.upper(), 0, i-1)
            rv_seq = ''
            i = len(seq)-1
            while i > -1:
                c = seq[i]
                if c.islower():
                    rv_seq += c
                else:
                    break
                i -= 1
            rv_primer = Primer(complement(rv_seq.upper()), i+1, len(seq)-1)
            if len(fw_seq) < 1 or len(rv_seq) < 1:
                raise RuntimeError("Error: Could not find lowercase sequence indicating amplicon primers.")
            primers[name] = [(fw_primer, rv_primer)]

    # apply one set of primers to all RNAs if needed
    if len(primers['__ALL_TARGETS__']) > 0:
        p = primers['__ALL_TARGETS__']
        for RNA in primers:
            primers[RNA] = [(deepcopy(x[0]),deepcopy(x[1])) for x in p]
    del primers['__ALL_TARGETS__']

    # locate primers by matching sequence
    for name, seq in iterate_fasta(fastas):
        seq = seq.upper()
        for n in range(len(primers[name])):
            if primers[name][n][0].left is None:
                # locate forward primer (match sequence)
                fp = primers[name][n][0].seq
                try:
                    left = seq.index(fp)
                    primers[name][n][0].left = left
                    primers[name][n][0].right = left + len(fp) - 1
                    try:
                        left2 = seq.index(fp, left+1)
                        raise RuntimeError("Multiple locations match forward primer in RNA '{}'. Explicit location input is not yet implemented.".format(name))
                    except ValueError:
                        pass
                except ValueError:
                    raise RuntimeError("Unable to locate forward primer in the sequence of RNA '{}'".format(name))
            if primers[name][n][1].left is None:
                # locate reverse primer (match reverse complement)
                rp = reverse_complement(primers[name][n][1].seq)
                try:
                    left = seq.index(rp)
                    primers[name][n][1].left = left
                    primers[name][n][1].right = left + len(rp) - 1
                    try:
                        left2 = seq.index(rp, left+1)
                        raise RuntimeError("Multiple locations match reverse primer in RNA '{}'. Explicit location input is not yet implemented.".format(name))
                    except ValueError:
                        pass
                except ValueError:
                    raise RuntimeError("Unable to locate reverse primer in the sequence of RNA '{}'".format(name))

    # check for targets with no matching primers
    for RNA in primers:
        if len(primers[RNA]) == 0:
            raise RuntimeError("Error: No amplicon primers provided for RNA '{}'".format(RNA))

    return primers


def write_primers(rna_name,
                  primer_pair_list,
                  filename,
                  n_pairs_filename):
    o = open(filename, 'w')
    o.write('>'+rna_name+'\n')
    for pair in primer_pair_list:
        o.write('{} {}\n'.format(pair[0].seq, pair[1].seq))
        o.write('{} {} {} {}\n'.format(pair[0].left,
                                       pair[0].right,
                                       pair[1].left,
                                       pair[1].right))
    o.close()

    o2 = open(n_pairs_filename, 'w')
    o2.write('{}'.format(len(primer_pair_list)))
    o2.close()


if __name__=="__main__":
    ap = ArgumentParser()
    ap.add_argument("--fastas", type=str, nargs='+', required=True)
    ap.add_argument("--target-names", type=str, nargs='+', required=True)
    ap.add_argument("--primer-files", type=str, nargs='+')
    ap.add_argument("--primers-in-sequence", action="store_true")
    ap.add_argument("--locations-out", type=str, nargs='+', required=True)
    ap.add_argument("--n-pairs-out", type=str, nargs='+', required=True)
    pa = ap.parse_args()
    if len(pa.target_names) != len(pa.locations_out):
        raise RuntimeError("Error: number of output files must match the number of RNA targets.")
    primers = locate_primer_pairs(pa.fastas,
                                  pa.target_names,
                                  pa.primer_files,
                                  pa.primers_in_sequence)
    for i in range(len(pa.target_names)):
        rna = pa.target_names[i]
        write_primers(rna,
                      primers[rna],
                      pa.locations_out[i],
                      pa.n_pairs_out[i])
