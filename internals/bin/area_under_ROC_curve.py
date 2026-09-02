#!/usr/bin/env python3
import sys
from argparse import ArgumentParser
import numpy as np

from sklearn.metrics import roc_curve, auc

def load_ct(filename):
    f = open(filename, "r")
    f.readline()
    paired = []
    for line in f:
        s = line.strip().split()
        if s[4] == '0':
            paired.append(0)
        else:
            paired.append(1)
    paired = np.array(paired, dtype=bool)
    unpaired = ~paired
    return unpaired


def load_map(filename):
    f = open(filename, "r")
    profile = []
    for line in f:
        s = line.strip().split()
        val = float(s[1])
        if val < -990:
            val  = np.nan
        profile.append(val)
    profile = np.array(profile)
    return profile


ap = ArgumentParser()
ap.add_argument("--ct", type=str)
ap.add_argument("--map", type=str)
ap.add_argument("--min-auc", type=float)
ap.add_argument("--name", type=str)
pa = ap.parse_args(sys.argv[1:])

unpaired = load_ct(pa.ct)
profile = load_map(pa.map)

good = np.isfinite(profile)
unpaired = unpaired[good]
profile = profile[good]

FPR, TPR, _ = roc_curve(unpaired, profile)
# TODO: optionally output tab-delimited table for comparison plots between versions

AUC = auc(FPR, TPR)
print("{} AUC: {:0.3f}".format(pa.name, AUC))
if (AUC < pa.min_auc) and (pa.min_auc-AUC > 0.0001):
    sys.exit(1)
else:
    sys.exit(0)
