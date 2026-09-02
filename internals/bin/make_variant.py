#!/usr/bin/env python3
"""
Inputs:
 - fasta file containing 1 sequence
   or fasta file containing multiple sequences + selected sequence name
 - file containing sequencing depths and variant counts

Given variants and depths, calc variant freqs, only keeping those above
a depth and frequency threshold
- Update sequence and output a single corrected fasta
- Warn about variants above depth threshold and within some range of
  intermediate frequencies (maybe between 10% and the upper freq threshold)
- In the future, do k-means clustering to group variants of similar frequency,
  and output a fasta with multiple sequences corresponding to major species

More rigorous approach would be to build covariance matrix and chain localized variants
into larger mutants

Most rigorous would be a full de-novo assembly (too much of an art to automate reliably)

"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os, sys, argparse
import numpy as np


fmt_msg = "Error: variant counts file appears misformatted."
fmt_err = RuntimeError(fmt_msg)

len_msg = "Error: number of lines in variant count file "
len_msg += "does not match number of lines in depth file."
len_err = RuntimeError(len_msg)


# TODO: move load_fasta(), load_depth() to their own file utility module, since
# they are used by both make_reactivity_profiles.py and make_variant.py
def load_fasta(fastaname, rna, convert_to_rna=True):
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
                s += " sequence was specified (use the --rna option)."
                raise RuntimeError(s)
            if name == rna:
                in_selected_rna = True
            seq_name = name
            continue
        if in_selected_rna:
            seq += line.strip().replace(' ','').replace('\t', '')
    if convert_to_rna:
        seq = seq.replace("T", "U")
    return seq_name, seq


class Variant:
    def __init__(self,
                 left,
                 right,
                 seq,
                 count):
        assert isinstance(left, int)
        assert isinstance(right, int)
        assert isinstance(seq, str)
        assert isinstance(count, int)
        self.left = left  # 0-based
        self.right = right  # 0-based
        self.seq = seq
        self.count = count

    def __str__(self):
        # convert to 1-based coords when printing
        return '({}-{}, "{}", {})'.format(self.left + 1,
                                          self.right + 1,
                                          self.seq,
                                          self.count)


class VariantFreq:
    def __init__(self,
                 left,
                 right,
                 seq,
                 freq):
        assert isinstance(left, int)
        assert isinstance(right, int)
        assert isinstance(seq, str)
        assert isinstance(freq, float)
        self.left = left  # 0-based
        self.right = right  # 0-based
        self.seq = seq
        self.freq = freq

    def __str__(self):
        # convert to 1-based coords when printing
        return '({}-{}, "{}", {:.3f})'.format(self.left + 1,
                                              self.right + 1,
                                              self.seq,
                                              self.freq)


class VariantFreqVerbose:
    def __init__(self,
                 left,
                 right,
                 ref_seq,
                 rep_seq,
                 freq):
        assert isinstance(left, int)
        assert isinstance(right, int)
        assert isinstance(ref_seq, str)
        assert isinstance(rep_seq, str)
        assert isinstance(freq, float)
        self.left = left  # 0-based
        self.right = right  # 0-based
        self.reference_seq = ref_seq
        self.replacement_seq = rep_seq
        self.freq = freq

    def __str__(self):
        # convert to 1-based coords when printing
        return '({}-{}, "{}", "{}", {:.3f})'.format(self.left + 1,
                                                    self.right + 1,
                                                    self.reference_seq,
                                                    self.replacement_seq,
                                                    self.freq)


def parse_variants(s):
    """
     generator yielding Variant objects parsed from string

    """

    open_brace = False
    internal = ""
    for c in s:
        if c == "(":
            if open_brace:
                raise fmt_err
            open_brace = True
            internal = ""
        elif c == ")":
            if not open_brace:
                raise fmt_err
            open_brace = False
            fields = [x.strip() for x in internal.split(',')]
            if len(fields) != 3:
                raise fmt_err
            bounds = fields[0].split('-')
            if len(bounds) != 2:
                raise fmt_err
            try:
                left, right = (int(x) for x in bounds)
                if right <= left:
                    raise fmt_err
            except ValueError:
                raise fmt_err
            seq = fields[1].replace('"', '')
            try:
                count = int(fields[2])
            except ValueError:
                raise fmt_err
            variant = Variant(left, right, seq, count)
            yield variant
        elif open_brace:
            internal += c


def load_variants(variant_filepath,
                  seq_len,
                  mindepth=None,
                  minfreq=None,
                  maxspan=None,
                  maxlength=None,
                  warningfreq=None):
    # NOTE: variants returned from this function will have counts replaced
    #       by frequencies (VariantFreq instead of Variant)

    # variants with apparent frequencies between warningfreq and minfreq
    # will be reported to the user through a warning, but will not affect
    # corrected output sequence

    assert isinstance(mindepth, int)
    assert isinstance(minfreq, float)
    assert isinstance(maxspan, int)
    assert isinstance(maxlength, int)
    assert isinstance(warningfreq, float)

    # load depths in first pass
    f = open(variant_filepath, "r")
    depths = []
    for line in f:
        try:
            depths.append(int(line.split()[0]))
        except ValueError:
            raise fmt_err
    depths = np.array(depths)
    f.seek(0)

    if len(depths) != seq_len:
        s = "Error: length of input sequence does not match length of input sequence variant file"
        raise RuntimeError(s)

    filtered_variants = []
    warning_variants = []

    i = 0
    for line in f:
        i += 1
        if i > len(depths):
            raise len_err
        for variant in parse_variants(line):
            # excised_seq = seq[variant.left+1:variant.right]

            # use min depth across bounds as
            # approximate number of spanning reads
            #d = min(depths[variant.left:variant.right + 1])

            # now just using min depth on left and right, since otherwise
            # intervening gaps will mess things up (no basecall in a gap, so doesn't
            # contribute to depth)
            d = min(depths[variant.left], depths[variant.right])

            freq = variant.count / d
            #print("{}, {}, {:0.2f}".format(variant, d, freq))
            if (variant.right - variant.left - 1 <= maxspan and
                        len(variant.seq) <= maxlength and
                        d >= mindepth and
                        freq >= minfreq):
                v = VariantFreq(variant.left, variant.right, variant.seq, freq)
                filtered_variants.append(v)
            elif (variant.right - variant.left - 1 <= maxspan and
                  len(variant.seq) <= maxlength and
                  d >= mindepth and
                  freq >= warningfreq and
                  freq < minfreq):
                v = VariantFreq(variant.left, variant.right, variant.seq, freq)
                warning_variants.append(v)

    if i < len(depths):
        raise len_err

    return filtered_variants, warning_variants


def merge_adjacent_variants(variants):
    """
    Merge adjacent variant sequences so that final correction
    doesn't fail during iterative mapping. Assumes variants are
    ordered by leftmost unchanged nuc from 5-prime to 3-prime.

    """
    if len(variants) < 1:
        return variants
    new_variants = [variants.pop(0)]
    for v in variants:
        if new_variants[-1].right - 1 == v.left:
            new_variants[-1].right = v.right
            new_variants[-1].seq += v.seq
        else:
            new_variants.append(v)
    return new_variants


def add_refseqs(seq, variants):
    """
    Add reference sequences over range of each mutation for better
    human-readable output.

    """
    verbose_variants = []
    for v in variants:
        verbose_variants.append(VariantFreqVerbose(v.left,
                                                   v.right,
                                                   seq[v.left+1:v.right],
                                                   v.seq,
                                                   v.freq))
    return verbose_variants


def correct_sequence(seq,
                     variants):
    """
    Given a sequence and a list of variants, generate a corrected
    sequence.
    """

    # assumes no adjacent variants (i.e. merge_adjacent_variants already run)

    # maintain an array of changed regions
    # to guard against conflicting variants
    # being silently merged
    changed_nucs = np.zeros(len(seq), dtype=bool)
    num_corrections = 0
    new_seq = str(seq)
    mapping = list(range(len(seq)))
    # print(seq)

    # - maintain mapping of current seq coords to original seq coords
    # - to convert original seq coords into current adjusted seq coords,
    #   just index() left and right
    for v in variants:
        for i in range(v.left + 1, v.right):
            if changed_nucs[i]:
                s = "Error: Unable to determine variant sequence. "
                s += "Conflicting (overlapping) variant sequences present."
                raise RuntimeError(s)
            changed_nucs[i] = True
        num_corrections += 1

        # - map to current adjusted sequence coords
        l = mapping.index(v.left)
        r = mapping.index(v.right)

        new_seq = new_seq[:l + 1] + v.seq + new_seq[r:]
        mapping = mapping[:l + 1] + [' '] * len(v.seq) + mapping[r:]

        # print(new_seq)

    new_seq = new_seq.replace("U", "T")

    return new_seq, num_corrections


def write_fasta(corrected_seq,
                seq_name,
                outpath,
                linewidth=80):
    f = open(outpath, "w")
    f.write(">" + seq_name + "\n")
    n = 0
    for i in range(len(corrected_seq)):
        if n > linewidth:
            n = 0
            f.write("\n")
        f.write(corrected_seq[i])
        n += 1


if __name__ == "__main__":

    # TODO: add full scanning mode for extremely large seqs?

    # TODO: move this to an actual testing module
    # test multiple adjacent variant merging
    # seq = "ATGCATGCATGCATGC"
    # # deletion of nuc 3
    # v1 = Variant(1, 3, '', 10)
    # # mismatch at nuc 4
    # v2 = Variant(2, 4, 'A', 10)
    # variants = [v1, v2]
    # print("variants")
    # for v in variants:
    #     print(str(v))
    # variants = merge_adjacent_variants(variants)
    # print("merged")
    # for v in variants:
    #     print(str(v))
    # new_seq, num_corrections = correct_sequence(seq, variants)
    # print("seq:")
    # print(seq)
    # print("corrected seq:")
    # print(new_seq)
    # exit()

    parser = argparse.ArgumentParser()

    h = "Fasta file containing sequence of RNA of interest."
    parser.add_argument("--fa", required=True, type=str, help=h)

    h = "Name of RNA of interest. Required argument if fasta"
    h += " file contains more than one sequence."
    parser.add_argument("--rna", type=str, help=h)

    h = "Counted sequence variant file."
    parser.add_argument("--variants", required=True, type=str, help=h)

    h = "Output fasta file with corrected variant sequence."
    parser.add_argument("--out", type=str, required=True, help=h)

    h = "Maximum consecutive nucleotides changed (original nucs)."
    parser.add_argument("--maxspan", type=int, default=1000, help=h)
    h = "Maximum consecutive nucleotides in any sequence correction (replacement nucs)."
    parser.add_argument("--maxlength", type=int, default=1000, help=h)
    h = "Minimum sequencing depth for any sequence correction."
    parser.add_argument("--mindepth", type=int, default=50, help=h)
    h = "Minimum frequency for any sequence correction."
    parser.add_argument("--minfreq", type=float, default=0.6, help=h)
    h = "Frequency threshold above which to warn the user (but not change output sequence)."
    parser.add_argument("--warningfreq", type=float, default=0.1, help=h)

    p = parser.parse_args(sys.argv[1:])

    seq_name, seq = load_fasta(p.fa, p.rna, convert_to_rna=False)

    filtered_variants, warning_variants = load_variants(p.variants,
                                                        len(seq),
                                                        mindepth=p.mindepth,
                                                        minfreq=p.minfreq,
                                                        maxspan=p.maxspan,
                                                        maxlength=p.maxlength,
                                                        warningfreq=p.warningfreq)

    # print("filtered variants:")
    # for v in filtered_variants:
    #    print(v)

    merged_variants = merge_adjacent_variants(filtered_variants)
    #merged_warning_variants = merge_adjacent_variants(warning_variants)

    corrected_seq, num_corrections = correct_sequence(seq,
                                                      merged_variants)

    verbose_merged_variants = add_refseqs(seq,
                                          merged_variants)
    verbose_warning_variants = add_refseqs(seq,
                                           warning_variants)

    write_fasta(corrected_seq,
                seq_name,
                p.out)

    msg = "{} correction{} made to sequence \"{}\"\n".format(num_corrections,
                                                             's' if num_corrections > 1 or num_corrections == 0 else '',
                                                             seq_name)
    if num_corrections > 0:
        msg += "Sequence changes:\n"
        msg += "(left-most 1-based unchanged nucleotide,\n"
        msg += "right-most 1-based unchanged nucleotide,\n"
        msg += "original sequence,\n"
        msg += "replacement sequence,\n"
        msg += "frequency)\n"
    sys.stdout.write(msg)
    for v in verbose_merged_variants:
        sys.stdout.write(" {}\n".format(v))
    sys.stdout.flush()

    if len(warning_variants) > 0:
        msg = "\n"
        msg += "WARNING: the following sequence changes have frequencies below {:.3f},\n".format(p.minfreq)
        msg += "but above {:.3f}:\n".format(p.warningfreq)
        msg += "(left-most 1-based unchanged nucleotide,\n"
        msg += "right-most 1-based unchanged nucleotide,\n"
        msg += "original sequence,\n"
        msg += "replacement sequence,\n"
        msg += "frequency)\n"
        sys.stdout.write(msg)
        for v in verbose_warning_variants:
            sys.stdout.write(" {}\n".format(v))
        msg = "These changes are not included the sequence written to {}\n".format(p.out)
        sys.stdout.write(msg)
        sys.stdout.flush()
