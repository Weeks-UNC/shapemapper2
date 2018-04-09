#!/usr/bin/env python

import sys, os, math
import unittest
from argparse import ArgumentParser

import numpy as np
from numpy import isnan, nan

""" Example scale factors file format:
A   0.474700
U   0.122603
G   0.272967
C   0.129731
"""

# derived from median unpaired background-subtracted reactivities
# for extracted E. coli ribosomes (both subunits)
global_scale_factors = {
    "5NIA":{
        'A':0.474700,
        'U':0.122603,
        'G':0.272967,
        'C':0.129731,
    },
    "1M7":{
        'A': 0.301910,
        'U': 0.269965,
        'G': 0.265406,
        'C': 0.162718,
    }
}


def load_scale_factors(filename):
    f = open(filename, "rU")
    nucs = "AUGC"
    scale_factors = {}
    for line in f:
        s = line.strip().split()
        if len(s) < 2:
            continue
        if s[0] not in nucs:
            raise RuntimeError("Scale factors file appears misformatted. First column should be nucs, A|U|G|C")
        try:
            val = float(s[1])
        except ValueError:
            raise RuntimeError("Scale factors file appears misformatted. Second column should be numeric scale factors")
        scale_factors[s[0]] = val
    if len(scale_factors) != 4:
        raise RuntimeError("Scale factors file does not contain all 4 nucs")
    return scale_factors


def load_map(filename):
    f = open(filename, "rU")
    nums = []
    seq = []
    profile = []
    stderrs = []
    for line in f:
        s = line.strip().split('\t')
        try:
            nums.append(s[0])
            seq.append(s[3].replace("T","U"))
            profile.append(float(s[1]))
            stderrs.append(float(s[2]))
            if profile[-1] < -990:
                profile[-1] = np.nan
            if stderrs[-1] < -990:
                stderrs[-1] = np.nan
        except (IndexError, ValueError) as e:
            raise RuntimeError("map file appears misformatted")
    return nums, seq, profile, stderrs


def write_map(nums, seq, profile, stderrs, filename):
    o = open(filename, "w")
    for i in range(len(nums)):
        val = "{:.6f}".format(profile[i])
        if not np.isfinite(profile[i]):
            val = "-999"
        err = "{:.6f}".format(stderrs[i])
        if not np.isfinite(stderrs[i]):
            err = "-999"
        o.write("{}\t{}\t{}\t{}\n".format(nums[i],
                                          val,
                                          err,
                                          seq[i]))
    print("Wrote map file to {}".format(filename))


def rescale(seq, vals, stderrs,
            scale_factors,
            how="down"):
    nucs = ["A","U","G","C"]
    sf = {nuc: 1/scale_factors[nuc] for nuc in nucs}
    if how == "down":
        max_f = max(sf.values())
        sf = {nuc:sf[nuc]/max_f for nuc in nucs}
    elif how == "up":
        min_f = min(sf.values())
        sf = {nuc:sf[nuc]/min_f for nuc in nucs}
    else:
        raise RuntimeError('Unrecognized "how" parameter for rescale_5NIA(). Options: "up", "down"')
    scaled_vals = np.full((len(vals),), np.nan)
    scaled_stderrs = np.full((len(vals),), np.nan)
    for i in range(len(seq)):
        try:
            f = sf[seq[i]]
        except KeyError:
            f = 1.0
        scaled_vals[i] = vals[i] * f
        scaled_stderrs[i] = stderrs[i] * f
    return scaled_vals, scaled_stderrs



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
    x = np.extract(np.isfinite(array), array)
    if x.shape[0] < 10:
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
    return norm_factor


def normalize_profile(profile, stderrs):
    norm_factor = find_boxplot_factor(np.array(profile))
    norm_profile = profile/norm_factor
    norm_stderrs = stderrs/norm_factor
    return norm_profile, norm_stderrs


class ScaleTest(unittest.TestCase):
    def runTest(self):
        seq = "AUGC#AUGC"
        profile = [1., 1., 1., 1., 5., 0.474700, 0.122603, 0.272967, 0.129731]
        scaled_up, _ = rescale(seq, profile, profile, global_scale_factors["5NIA"], how="up")
        scaled_down, _ = rescale(seq, profile, profile, global_scale_factors["5NIA"], how="down")
        #print("seq: {}".format(seq))
        #print("profile: {}".format(profile))
        #print("scaled up: {}".format(scaled_up))
        #print("scaled down: {}".format(scaled_down))
        self.assertAlmostEqual(scaled_up[0], 1.)
        self.assertAlmostEqual(scaled_down[1], 1.)


class EndToEndTests(unittest.TestCase):
    def runTest(self):
        from subprocess import check_call
        from os.path import abspath
        from os import remove
        THIS_SCRIPT = abspath(__file__)

        # write dummy map file for input
        open("__dummy.map__", "w").write("""1\t-999\t0\t5NIA
1\t0.474700\t0.1\tA
2\t0.122603\t0.1\tU
3\t0.272967\t0.1\tG
4\t0.129731\t0.1\tC
5\t-999\t0\t1M7
6\t0.301910\t0.1\tA
7\t0.269965\t0.1\tU
8\t0.265406\t0.1\tG
9\t0.162718\t0.1\tC
10\t-999\t0\tN
11\t1.0000\t0.1\tA
12\t1.0000\t0.1\tU
13\t1.0000\t0.1\tG
14\t1.0000\t0.1\tC
15\t-999\t0\tA""")

        # write dummy scale factors to file for optional input
        open("__dummy_scale_factors.txt__", "w").write("""A   0.474700
U   0.122603
G   0.272967
C   0.129731
""")

        cmd = "python {} --reagent 1M7 __dummy.map__ __dummy_scaled_1M7.map__".format(THIS_SCRIPT)
        print("Running cmd:")
        print(cmd)
        check_call(cmd, shell=True)

        cmd = "python {} --no-renorm __dummy.map__ __dummy_scaled_5NIA_no_renorm.map__".format(THIS_SCRIPT)
        print("Running cmd:")
        print(cmd)
        check_call(cmd, shell=True)

        cmd = "python {} --scale-factors __dummy_scale_factors.txt__ __dummy.map__ __dummy_scaled_5NIA_external.map__".format(THIS_SCRIPT)
        print("Running cmd:")
        print(cmd)
        check_call(cmd, shell=True)

        # remove dummy files
        remove("__dummy.map__")
        remove("__dummy_scaled_1M7.map__")
        remove("__dummy_scaled_5NIA_no_renorm.map__")
        remove("__dummy_scaled_5NIA_external.map__")
        remove("__dummy_scale_factors.txt__")


def run_tests():
    suite = unittest.TestSuite()
    suite.addTest(ScaleTest())
    suite.addTest(EndToEndTests())
    runner = unittest.TextTestRunner()
    runner.run(suite)


def check_test_arg(args):
    tp = ArgumentParser(add_help=False)
    tp.add_argument("--test", action="store_true")
    pta, rest = tp.parse_known_args(sys.argv[1:])
    if pta.test:
        run_tests()
        sys.exit()



if __name__=="__main__":
    check_test_arg(sys.argv[1:])

    ap = ArgumentParser()
    ap.add_argument("input", type=str, help="Input .map file")
    ap.add_argument("output", type=str, help="Output .map file")
    ap.add_argument("--reagent", type=str, default="5NIA", help="5NIA | 1M7")
    ap.add_argument("--scale-factors", type=str, default=None, help="Provide scale factors in a separate file. See comment at the top of this script for example format.")
    ap.add_argument("--no-renorm", action="store_true", help="Do not renormalize profile after scaling.")

    if len(sys.argv[1:]) == 0:
        ap.print_help()
        ap.exit()

    pa = ap.parse_args(sys.argv[1:])

    nums, seq, profile, stderrs = load_map(pa.input)
    if pa.scale_factors is not None:
        scale_factors = load_scale_factors(pa.scale_factors)
    else:
        try:
            scale_factors = global_scale_factors[pa.reagent]
        except KeyError:
            raise RuntimeError("Precomputed scale factors not available for reagent \"{}\"".format(pa.reagent))

    profile, stderrs = rescale(seq, profile, stderrs,
                               scale_factors)

    if not pa.no_renorm:
        profile, stderrs = normalize_profile(profile, stderrs)

    write_map(nums, seq, profile, stderrs, pa.output)
