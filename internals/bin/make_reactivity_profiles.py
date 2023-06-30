#!/usr/bin/env python3
"""
Given:
 fasta file containing 1 sequence
 or fasta file containing multiple sequences + selected sequence name
 counts + depth files for 1-3 samples

 Calculate:
 simple reactivity profile
 reactivity profile, excluding high-background positions and high-BG positions

 Output:
 summary tab-delimited file

"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys, os
import argparse

import numpy as np
from numpy import sqrt, square, isnan, nan

np.seterr("ignore") # disable divide-by-zero warnings

# counts and depths should be dicts of numpy arrays of int
# where dict keys are: modified, untreated, denatured
def calc_profile(counts, depths):
    # adjust read depth for regions where multinuc mutations prevent
    # resolution of modifications within their length

    rate_rx = None
    rate_bg = None
    rate_dc = None

    # single sample (don't use controls)
    if counts["denatured"] is None and counts["untreated"] is None:
        rate_rx = counts["modified"] / depths["modified"]  # Note: requires python3 (otherwise need np.true_divide())
        profile = np.copy(rate_rx)
        stderrs = sqrt(rate_rx) / sqrt(depths["modified"])

    # two samples (use untreated control)
    elif counts["denatured"] is None:
        rate_rx = counts["modified"] / depths["modified"]
        stderr_rx = sqrt(rate_rx) / sqrt(depths["modified"])
        rate_bg = counts["untreated"] / depths["untreated"]
        stderr_bg = sqrt(rate_bg) / sqrt(depths["untreated"])
        profile = rate_rx - rate_bg
        # sqrt of sum of the squares of the individual stderrs
        stderrs = sqrt(square(stderr_rx) + square(stderr_bg))

    # three samples (use both controls)
    else:
        rate_rx = counts["modified"] / depths["modified"]
        stderr_rx = sqrt(rate_rx) / sqrt(depths["modified"])
        rate_bg = counts["untreated"] / depths["untreated"]
        stderr_bg = sqrt(rate_bg) / sqrt(depths["untreated"])
        rate_dc = counts["denatured"] / depths["denatured"]
        stderr_dc = sqrt(rate_dc) / sqrt(depths["denatured"])
        diff = rate_rx - rate_bg
        profile = diff / rate_dc
        stderrs = sqrt(square(stderr_rx / rate_dc) +
                       square(stderr_bg / rate_dc) +
                       square(stderr_dc * diff / square(rate_dc)))

    rates = {"modified": rate_rx,
             "untreated": rate_bg,
             "denatured": rate_dc}

    # replace any non-finite values with NANs (otherwise sometimes get some "-inf", etc)
    for i in range(len(profile)):
        if not np.isfinite(profile[i]):
            profile[i] = nan
        if not np.isfinite(stderrs[i]):
            stderrs[i] = nan

    return profile, stderrs, rates

# exclude nucs with high background mutation rates
# exclude nucs with low read coverage
# exclude nucs with N in the target seq
#
# depths should be a dict of numpy arrays of int
# where dict keys are: modified, untreated, denatured
def filter_profile(profile, stderrs, depths, rates,
                   max_background=0.05, min_depth=1000,
                   random_primer_len=None):

    filtered_profile = np.copy(profile)
    filtered_stderr = np.copy(stderrs)

    for i in range(len(profile)):
        # exclude left-most nuc (since mutations on the end of sequence aren't reliable
        if i == 0:
            good_nuc = False
        if random_primer_len is not None and random_primer_len>0:
            if len(profile)-i <= random_primer_len+1:
                good_nuc = False
        good_nuc = True
        if seq[i].upper() not in ['A','U','G','C']:
            good_nuc = False
        if seq[i].islower():
            good_nuc = False
        if rates["untreated"] is not None:
            if rates["untreated"][i] > max_background:
                good_nuc = False
        for k in ["modified", "untreated", "denatured"]:
            if depths[k] is not None:
                if depths[k][i] < min_depth:
                    good_nuc = False

        if not isnan(profile[i]) and good_nuc:
            filtered_profile[i] = profile[i]
            filtered_stderr[i] = stderrs[i]
        else:
            filtered_profile[i] = nan
            filtered_stderr[i] = nan
    return filtered_profile, filtered_stderr


def write_tab_file(seq,
                   counts, read_depths, effective_depths, rates,
                   profile, stderrs,
                   filtered_profile, filtered_stderrs,
                   addtl_data, addtl_column_names,
                   out,
                   decimal_places=6):
    n = "{{:.{}f}}".format(decimal_places)

    headers_types = [
        ('Nucleotide','{}'),
        ('Sequence','{}'),

        ('Modified_mutations','{}'),
        ('Modified_read_depth', '{}'),
        ('Modified_effective_depth','{}'),
        ('Modified_rate',n),
    ]
    for i, name in enumerate(addtl_column_names):
        headers_types += [('Modified_'+name, '{}')]

    headers_types += [
        ('Untreated_mutations','{}'),
        ('Untreated_read_depth', '{}'),
        ('Untreated_effective_depth','{}'),
        ('Untreated_rate',n),
    ]
    for i, name in enumerate(addtl_column_names):
        headers_types += [('Untreated_'+name, '{}')]

    headers_types += [
        ('Denatured_mutations','{}'),
        ('Denatured_read_depth', '{}'),
        ('Denatured_effective_depth','{}'),
        ('Denatured_rate',n),
    ]
    for i, name in enumerate(addtl_column_names):
        headers_types += [('Denatured_'+name, '{}')]

    headers_types += [
        ('Reactivity_profile',n),
        ('Std_err',n),

        ('HQ_profile',n),
        ('HQ_stderr',n),
    ]
    headers, types = zip(*headers_types)

    f = open(out, "w")
    f.write("\t".join(headers)+"\n")

    for i in range(len(seq)):
        line = "\t".join(types)
        fields = [
            i+1,
            seq[i],
        ]
        for sample in ("modified","untreated","denatured"):
            if rates[sample] is not None:
                fields += [counts[sample][i],
                           read_depths[sample][i],
                           effective_depths[sample][i],
                           rates[sample][i]]
            else:
                fields += [0, 0, 0, nan]
            for n in range(len(addtl_column_names)):
                if addtl_data[sample] is not None:
                    fields += [addtl_data[sample][n][i]]
                else:
                    fields += [0]

        fields += [profile[i],
                   stderrs[i],

                   filtered_profile[i],
                   filtered_stderrs[i]
                   ]
        line = line.format(*fields)
        f.write(line)
        if i < len(seq)-1:
            f.write("\n")

def load_fasta(fastaname, rna, convert_to_rna=True):
    f = open(fastaname, "rU")
    seq = ""
    rna_count = 0
    in_selected_rna = False
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
            continue
        if in_selected_rna:
            seq += line.strip().replace(' ','').replace('\t', '')
    if convert_to_rna:
        seq = seq.replace("T", "U").replace('t','u')
    return seq


def load_counts(filename):
    if filename is None:
        return None

    f = open(filename, "rU")

    # do one pass to determine array length
    # TODO: might actually be faster to just resize array in memory and read in one pass
    length = 0
    f.readline()  # skip header
    for line in f:
        length += 1
    f.seek(0)

    if length == 0:
        s = "Error: mutation counts file "
        s += "\"" + filename + "\""
        s += " contains no data."
        raise RuntimeError(s)

    headers = f.readline().strip().split('\t')

    mismatch_headers = [
        "AT", "AG", "AC", "TA", "TG", "TC", "GA", "GT", "GC", "CA", "CT", "CG",
        "A_multinuc_mismatch", "C_multinuc_mismatch", "G_multinuc_mismatch", "T_multinuc_mismatch", "N_multinuc_mismatch"
    ]

    deletion_headers = [
        "A-", "T-", "G-", "C-",
        "multinuc_deletion",
        "complex_deletion",
    ]

    insert_headers = [
        "-A", "-T", "-G", "-C", "-N",
        "multinuc_insertion",
        "complex_insertion",
    ]

    read_depth_header = "read_depth"
    effective_depth_header = "effective_depth"

    expected_headers = list(mismatch_headers)
    expected_headers += deletion_headers
    expected_headers += insert_headers

    required_headers = list(expected_headers)

    # also load ambig columns if present
    expected_headers += [h+"_ambig" for h in mismatch_headers]
    expected_headers += [h+"_ambig" for h in deletion_headers]
    expected_headers += [h+"_ambig" for h in insert_headers]

    missing_headers = []
    for h in required_headers + \
            [read_depth_header] + \
            [effective_depth_header]:
        if h not in headers:
            missing_headers += [h]
    if len(missing_headers)>0:
        s = "Mutation counts file "
        s += "\"" + filename + "\""
        s += " is missing expected header(s): " + ', '.join(missing_headers)
        raise RuntimeError(s)

    # mapped depth input/output
    '''
    "off_target_mapped_depth"
    "low_mapq_mapped_depth"
    "primer_pair_n_mapped_depth" - used with --amplicon
    "mapped_depth" - used if no --amplicon
    '''
    addtl_columns = []
    addtl_column_names = []
    for i, h in enumerate(headers):
        if (h in ["off_target_mapped_depth",
                  "low_mapq_mapped_depth",
                  "mapped_depth"] or
            (h.startswith("primer_pair_") and
             h.endswith("_mapped_depth"))):
            addtl_column_names.append(h)
            addtl_columns.append(i)


    selected_columns = []
    selected_column_names = []
    for h in expected_headers:
        if h in headers:
            selected_columns.append(headers.index(h))
            selected_column_names.append(h)
    read_depth_column = headers.index(read_depth_header)
    effective_depth_column = headers.index(effective_depth_header)

    # now load counts and depths
    counts = np.zeros(length, dtype=np.int)
    read_depths = np.zeros(length, dtype=np.int)
    effective_depths = np.zeros(length, dtype=np.int)
    addtl_data = []
    for i in range(len(addtl_column_names)):
        addtl_data.append(np.zeros(length, dtype=np.int))
    try:
        for i in range(length):
            line = f.readline().strip()
            s = line.split('\t')
            for n in selected_columns:
                counts[i] += int(s[n])
            read_depths[i] += int(s[read_depth_column])
            effective_depths[i] += int(s[effective_depth_column])
            for ni, n in enumerate(addtl_columns):
                addtl_data[ni][i] += int(s[n])

    except IndexError:
        s = "Mutation counts file "
        s += "\"" + filename + "\""
        s += " is misformatted (line contains too few columns)."
        raise RuntimeError(s)
    except ValueError:
        s = "Mutation counts file "
        s += "\"" + filename + "\""
        s += " is misformatted (contains non-numeric characters)."

    print("loaded counts from these columns in mutation counts file: {}".format(selected_column_names))
    print("loaded read depths from column {} in mutation counts file".format(read_depth_header))
    print("loaded effective read depths from column {} in mutation counts file".format(effective_depth_header))

    return counts, read_depths, effective_depths, addtl_data, addtl_column_names


if __name__ == "__main__":
    # fasta file containing 1 sequence
    # or fasta file containing multiple sequences + selected sequence name
    # counts + depth files for 1-3 samples

    # TODO: add some sort of scanning mode for extremely large seqs?

    parser = argparse.ArgumentParser()

    h = "Fasta file containing sequence of RNA of interest."
    h += " DNA sequences will be converted to RNA automatically."
    h += " Lowercase regions will be excluded from profile calculation."
    parser.add_argument("--fa", required=True, type=str, help=h)

    h = "Name of RNA of interest. Required argument if fasta"
    h += " file contains more than one sequence."
    parser.add_argument("--rna", type=str, help=h)

    h = "Counted mutations files in the following order:"
    h += " modified, untreated control, denatured control."
    h += " If the denatured control is omitted, it will not"
    h += " contribute to reactivity profile calculation. If"
    h += " only modified counts are provided, neither control"
    h += " will contribute to reactivity profile calculation."
    parser.add_argument("--counts", required=True, type=str, nargs="*", help=h)

    h = "Output tab-delimited file."
    parser.add_argument("--out", type=str, required=True, help=h)

    parser.add_argument("--mindepth", type=int, default=None)
    parser.add_argument("--maxbg", type=float, default=None)
    parser.add_argument("--random-primer-len", type=int, default=None)
    parser.add_argument("--dms", action='store_true', help='Use DMS normalization')
    

    p = parser.parse_args(sys.argv[1:])

    print("argv: "+str(sys.argv))

    if len(p.counts) > 3:
        raise ValueError("Error: too many mutation counts files provided (3 max).")

    seq = load_fasta(p.fa, p.rna)
    if p.rna is not None:
        print("loaded sequence {} from file {}".format(p.rna, p.fa))
    else:
        print("loaded sequence from file {}".format(p.fa))

    samples = ["modified", "untreated", "denatured"]
    counts = {}
    read_depths = {}
    effective_depths = {}
    addtl_data = {}
    addtl_column_names = []

    for i in range(len(samples)):
        try:
            counts[samples[i]], \
            read_depths[samples[i]], \
            effective_depths[samples[i]], \
            addtl_data[samples[i]], \
            addtl_column_names = load_counts(p.counts[i])
            print("loaded counts and depths for {} sample from file {}".format(samples[i], p.counts[i]))
        except IndexError:
            counts[samples[i]] = None
            read_depths[samples[i]] = None
            effective_depths[samples[i]] = None
            addtl_data[samples[i]] = None
        
    # check that seq length matches all mutation count data length and depth length
    lengths = []
    lengths.append(len(seq))
    for k in samples:
        if counts[k] is not None:
            lengths.append(counts[k].shape[0])
        if read_depths[k] is not None:
            lengths.append(read_depths[k].shape[0])
        if effective_depths[k] is not None:
            lengths.append(effective_depths[k].shape[0])
    if len(set(lengths)) > 1:
        s = "Error: input data lengths do not all match."
        raise RuntimeError(s)

    profile, stderrs, rates = calc_profile(counts, effective_depths)

    args = profile, stderrs, effective_depths, rates
    kwargs = {}
    if p.mindepth is not None:
        kwargs["min_depth"] = p.mindepth
    if p.maxbg is not None:
        kwargs["max_background"] = p.maxbg
    if p.random_primer_len is not None:
        kwargs["random_primer_len"] = p.random_primer_len
    filtered_profile, filtered_stderrs = filter_profile(*args, **kwargs)

    # create output directory if needed
    d = os.path.split(p.out)[0]
    if len(d)>0:
        os.makedirs(d, exist_ok=True)

    write_tab_file(seq,
                   counts, read_depths, effective_depths, rates,
                   profile, stderrs,
                   filtered_profile, filtered_stderrs,
                   addtl_data, addtl_column_names,
                   p.out)
    print("wrote profile to {}".format(p.out))
