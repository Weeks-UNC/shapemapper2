#!/usr/bin/env python3
"""
 Given:
 normalization targets:
  a single tab-delimited reactivity profile summary file
  -or-
  a list of such files
 profiles to scale:
  no files (in which case the input files will be normalized)
  -or-
  a tab file or files

 Calculate:
 a normalization (scaling) factor using all normalization targets

 Output:
 overwrite profiles with new normalized data column
 -or-
 if --normout or --scaleout files are specified, output files will be
 produced instead of overwriting input files
"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

# TODO: add some sort of scanning mode for extremely large seqs?
# TODO: support .shape and .map files
# FIXME: in overwrite mode write to a temporary file, then replace input file
#        (in case script crashes mid-write, file will not be lost)

import sys, os, math
import argparse

import numpy as np
from numpy import isnan, nan

#Following 3 functions modified from Gregg Rice's boxplot normalization script
#
# 1.find the scaling factor by ranking  all of the shape values
# 2. take 1.5*abs(Q1-Q3) as a cutoff
# 3. remove either the top 10% of the RNA or the positions above this cutoff, whichever is smaller
# 4. Average the next 10% from the original length of the RNA --> this is the scaling factor
def calc_quartile(x, q, qtype=7):
    # source: http://adorio-research.org/wordpress/?p=125
    # x = array, q = quartile (in % as a decimal)
    y = np.copy(x)
    n = len(y)
    abcd = [(0, 0, 1, 0),  # inverse empirical distrib.function., R type 1
            (0.5, 0, 1, 0),  # similar to type 1, averaged, R type 2
            (0.5, 0, 0, 0),  # nearest order statistic,(SAS) R type 3
            (0, 0, 0, 1),  # California linear interpolation, R type 4
            (0.5, 0, 0, 1),  # hydrologists method, R type 5
            (0, 1, 0, 1),  # mean-based estimate(Weibull method), (SPSS,Minitab), type 6
            (1, -1, 0, 1),  # mode-based method,(S, S-Plus), R type 7
            (1.0 / 3, 1.0 / 3, 0, 1),  # median-unbiased ,  R type 8
            (3 / 8.0, 0.25, 0, 1)  # normal-unbiased, R type 9.
            ]
    a, b, c, d = abcd[qtype - 1]
    g, j = math.modf(a + (n + b) * q - 1)
    if j < 0:
        return x[0]
    elif j >= n:
        return x[n - 1]
    j = int(math.floor(j))
    if g == 0:
        return x[j]
    else:
        return y[j] + (y[j + 1] - y[j]) * (c + d * g)


class NormError(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)


def find_boxplot_factor(array):
    x, o, a = [], [], 0
    # Following deprecated line is behavior that normalization and
    # structure modeling were optimized with, but this behavior
    # is probably not ideal. For RNAs with regions of poor sequencing
    # depth, treating those regions as unreactive could bias the
    # normalized reactivities. This is especially important for
    # larger RNAs.
    # x = np.fromiter((n if not isnan(n) else 0 for n in array))
    x = array[np.where(np.isfinite(array))]
    if x.shape[0] < 20:
        s = "Error: sequence contains too few nucleotides"
        s += " with quality reactivity information for"
        s += " effective normalization factor calculation."
        raise NormError(s)
    else:
        x.sort()
        ten_pct = len(x) // 10
        five_pct = len(x) // 20
        # calculate the interquartile range *1.5
        q_limit = 1.5 * abs(calc_quartile(x, 0.25) - calc_quartile(x, 0.75))
        ten_limit = x[x.shape[0] - 1 - ten_pct]
        five_limit = x[x.shape[0] - 1 - five_pct]
        # choose the cutoff that eliminates the fewest points
        limit = max(q_limit, ten_limit)
        if len(x) < 100:
            limit = max(q_limit, five_limit)
        # make new list without the outliers
        for i in range(len(x)):
            if x[i] < limit:
                o.append(x[i])
        # avg next ten percent
        try:
            for i in range(-ten_pct, 0):
                a = o[i] + a
            norm_factor = a / ten_pct
        except IndexError:
            raise NormError("Unable to calculate a normalization factor.")

    norm_factor_dict = {}
    for n in ('A','C','U', 'G'):
        norm_factor_dict[n] =  norm_factor

    return norm_factor_dict


def find_pernt_factor(sequence, profile):
    """This function performs per-nt reactivity normalization that is slightly different than boxplot normalization"""
    
    norm_factor_dict = {}

    for n in ('A','C','U','G'):

        mask = (sequence == n) & np.isfinite(profile)
        
        x = profile[mask]
        if x.shape[0] < 20:
            s = "Error: sequence contains too few "+n+" nucleotides (={})".format(x.shape[0])
            s += " with quality reactivity information for"
            s += " effective per-nt normalization factor calculation"
            raise NormError(s)

        else:
            bnds = np.percentile(x, [90., 95.])
            pmask = (x >= bnds[0]) & (x<bnds[1])
            normset = x[pmask]
            
            # compute norm standard way
            n1 = np.mean(normset)
            
            try:
                # compute the norm only considereing reactive nts
                n2 = np.percentile(x[x>0.001], 75)
            except IndexError:
                n2 = 0

            norm_factor_dict[n] = max(n1, n2)       
    

    # perform quality checks; if signals are too low, set to nan
    for n in norm_factor_dict:
        if norm_factor_dict[n] < 0.002:
            print('Signal for {} = {:.3f} is below acceptable level. All {} nts ignored'.format(n, norm_factor_dict[n], n))
            norm_factor_dict[n] = np.nan

    return norm_factor_dict   
    


def calc_norm_factor(sequences, profiles, pernt=False):
    
    combined_sequence = np.hstack(sequences)
    combined_profile = np.hstack(profiles)

    if pernt:
        return find_pernt_factor(combined_sequence, combined_profile)
    else:
        return find_boxplot_factor(combined_profile)


def normalize_profile(sequence, profile, stderrs, norm_factor):
    
    norm_profile = np.empty(profile.shape)
    norm_stderrs = np.empty(profile.shape)
    norm_profile[:] = np.nan
    norm_stderrs[:] = np.nan
    
    for i in range(profile.shape[0]):
        n = sequence[i]
        
        if n in norm_factor and np.isfinite(norm_factor[n]):
            norm_profile[i] = profile[i]/norm_factor[n]
            norm_stderrs[i] = stderrs[i]/norm_factor[n]
        else:
            norm_profile[i] = np.nan
            norm_stderrs[i] = np.nan

    return norm_profile, norm_stderrs



def load_profile(filename):
    
    f = open(filename)

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
    expected_fields = ['Sequence', 'HQ_profile','HQ_stderr']
    for x in expected_fields:
        if x not in headers:
            raise RuntimeError("File \""+filename+"\" does not contain the expected columns.")
    seq_col = headers.index(expected_fields[0])
    profile_col = headers.index(expected_fields[1])
    stderr_col = headers.index(expected_fields[2])
    
    sequence = []
    profile = np.empty(length)
    stderrs = np.empty(length)

    for i in range(length):
        line = f.readline()
        s = line.strip().split('\t')
        sequence.append(s[seq_col])
        profile[i] = float(s[profile_col])
        stderrs[i] = float(s[stderr_col])

    return np.array(sequence), profile, stderrs

# overwrite existing norm profile and stderr columns if present
# otherwise, add new columns. If an output file is given, write
# there instead of overwriting input file.
def write_norm_columns(profile,
                       stderrs,
                       filename,
                       outname,
                       decimal_places=6):
    n = "{{:.{}f}}".format(decimal_places)

    f = open(filename)
    lines = f.readlines()
    f.close()

    norm_prof_header = "Norm_profile"
    norm_stderr_header = "Norm_stderr"

    headers = lines[0].strip().split('\t')
    try:
        norm_profile_col = headers.index(norm_prof_header)
        norm_stderr_col = headers.index(norm_stderr_header)
    except ValueError:
        norm_profile_col = None
        norm_stderr_col = None

    if outname is None:
        f = open(filename, "w")
    else:
        # create path to output file if needed
        folder, filename = os.path.split(outname)
        if len(folder) > 0:
            os.makedirs(folder, exist_ok=True)
        f = open(outname, "w")
    
    lines.pop(0)
    if not ("Norm_profile" in headers or "Norm_stderr" in headers):
        headers.extend(["Norm_profile", "Norm_stderr"])
    f.write("\t".join(headers)+"\n")


    for i in range(len(lines)):
        line = lines[i]
        s = line.strip().split('\t')
        s1 = n.format(profile[i])
        if norm_profile_col is None:
            s.append(s1)
        else:
            s[norm_profile_col] = s1
        s2 = n.format(stderrs[i])
        if norm_stderr_col is None:
            s.append(s2)
        else:
            s[norm_stderr_col] = s2
        f.write("\t".join(s))
        if i < len(lines)-1:
            f.write("\n")


def dup(filename, outname):
    """
    Duplicate a text file. Just used to workaround the case
    when normalization can't be completed, but an output
    file is still expected.
    """
    f = open(filename)
    lines = f.readlines()
    f.close()
    o = open(outname, "w")
    o.writelines(lines)
    o.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    h = "List of tab-delimited file(s) from which to calculate "
    h += "reactivity profile normalization factor. Profiles will"
    h += " each be scaled using this factor and two additional"
    h += " columns will be added to these files, unless a list of"
    h += " output files is provided using --normout."
    parser.add_argument("--tonorm", type=str, nargs="*", required=True, help=h)

    h = "List of files to output. If this argument is specified,"
    h += " files given with --tonorm will not be overwritten."
    parser.add_argument("--normout", type=str, nargs="*", help=h)

    h = "List of additional tab-delimited file(s) to be scaled by the"
    h += "calculated normalization factor."
    parser.add_argument("--toscale", type=str, nargs="*", help=h)

    h = "List of files to output. If this argument is specified,"
    h += " files given with --toscale will not be overwritten."
    parser.add_argument("--scaleout", type=str, nargs="*", help=h)

    h = "Warn if normalization fails, instead of exiting with error."
    parser.add_argument("--warn-on-error", action="store_true", default=False, help=h)

    h = "Normalize data using per-nt DMS normalization scheme."
    parser.add_argument("--dms", action="store_true", default=False, help=h)

    p = parser.parse_args(sys.argv[1:])

    if p.tonorm is not None and p.normout is not None:
        if len(p.tonorm) != len(p.normout):
            msg = "Number of output files given with --normout must "
            msg += "match the number of input files given with --tonorm."
            raise RuntimeError(msg)

    if p.toscale is not None and p.scaleout is not None:
        if len(p.toscale) != len(p.scaleout):
            msg = "Number of output files given with --scaleout must "
            msg += "match the number of input files given with --toscale."
            raise RuntimeError(msg)


    if p.toscale is None:
        p.toscale = []

    if p.scaleout is None:
        p.scaleout = list(p.toscale)

    for i in range(len(p.tonorm)):
        f = p.tonorm[i]
        try:
            o = p.normout[i]
        except (AttributeError, TypeError) as e:
            # if no output files explicitly specified, overwrite columns in the input file
            o = None
        if f not in p.toscale:
            p.toscale.append(f)
            p.scaleout.append(o)
    
    sequence_list = []
    profile_list = []
    stderr_list = []

    for name in p.tonorm:
        sequence, profile, stderr = load_profile(name)
        print("loaded profile and stderrs from "+name)
        sequence_list.append(sequence)
        profile_list.append(profile)
        stderr_list.append(stderr)

    try:
        factors = calc_norm_factor(sequence_list, profile_list, pernt=p.dms)
        print("calculated normalization factor: {}".format(factors))
    
    except NormError as e:
        if p.warn_on_error:
            print(e)
            factors = None
        else:
            raise e

    for i in range(len(p.toscale)):
        input_filename = p.toscale[i]
        output_filename = p.scaleout[i]
        # FIXME: This ends up reading each file twice
        sequence, profile, stderr = load_profile(input_filename)
        if factors is not None:
            profile, stderr = normalize_profile(sequence, profile, stderr, factors)
            write_norm_columns(profile, stderr, input_filename, output_filename)
            #if output_filename is None:
            #    print("updated {} with normalized data columns".format(input_filename))
            #else:
            #    print("wrote normalized data from file {} to output file {}".format(input_filename, output_filename))
        else:
            # just copy data to new files, don't create new columns
            if output_filename is not None:
                dup(input_filename, output_filename)
