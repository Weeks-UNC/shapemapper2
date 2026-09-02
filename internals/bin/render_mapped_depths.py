#!/usr/bin/env python3
"""
Plot simple read coverage for each sample, broken down by amplicon
primer pairs if present, and showing reads excluded for low mapping
quality or off-target location.

"""


# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys, os, re
from math import floor, log10
from argparse import ArgumentParser as AP
from collections import OrderedDict
from cycler import cycler
import numpy as np

import matplotlib as mp
mp.use('Agg')
mp.rcParams["font.sans-serif"].insert(0,"Arial")
p = {"font.family": "sans-serif",
     "pdf.fonttype": 42, # use TrueType fonts when exporting PDFs
                         # (embeds most fonts - this is especially
                         #  useful when opening in Adobe Illustrator)

     'font.size': 10.0,
     'figure.titlesize': 'large',
     'figure.dpi': 100.0,
     'savefig.dpi': 100.0,     
     'legend.fancybox': True,
     'legend.shadow': False,
     'legend.edgecolor': '0.8',
     'lines.linewidth': 1.5,

     'xtick.direction': 'out',
     'ytick.direction': 'out',
     'legend.fontsize': 8,
     'grid.color': ".8",
     'grid.linestyle': '-',
     'grid.linewidth': 1,
     'axes.edgecolor': '.0',
     'axes.labelcolor': '.0',
     'xtick.color': '.0',
     'xtick.major.size': 6.0,
     'ytick.color': '.0',
     'ytick.major.size': 6.0,
     'xtick.minor.size': 2.0,
     'ytick.minor.size': 2.0,
     # color cycle from matplotlib 2
     'axes.prop_cycle': cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                         '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                         '#bcbd22', '#17becf']),
}
mp.rcParams.update(p)

mp.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


categories = ["included","excluded"]
samples = ["Modified","Untreated","Denatured"]


def metric_abbreviate(num):
    if isinstance(num, str):
        suffixes = {3:'k',
                    6:'M',
                    9:"G"}
        # replace trailing zeros with metric abbreviation
        zero_count = len(num)-len(num.rstrip('0'))
        suffix = ''
        new_string = str(num)
        for num_zeros in sorted(suffixes.keys()):
            if num_zeros <= zero_count:
                suffix = suffixes[num_zeros]
                new_string = num[:-num_zeros]
        new_string = new_string+suffix
    else:
        suffixes = ['','k','M','G',]
        num = float(num)
        d = max(0, min(len(suffixes)-1,
                       int(floor(0 if num == 0 else log10(abs(num))/3))))
        new_string = '{:.0f}{}'.format(num / 10**(3*d), suffixes[d])
    return new_string


def commify(s):
    o = ""
    for i in range(-1,-len(s)-1,-1):
        if -i%3==0:
            o = ','+s[i]+o
        else:
            o = s[i]+o
    o = o.lstrip(',')
    return o


class Primer:
    def __init__(self):
        self.index = 0
        self.seq = ""
        self.strand = '+'
        self.left = 0
        self.right = 0

    def __str__(self):
        return "{} {} {} {}".format(self.seq, self.strand, self.left, self.right)


class PrimerPair:
    def __init__(self):
        self.fw = Primer()
        self.fw.strand = '+'
        self.rv = Primer()
        self.rv.strand = '-'

    def __str__(self):
        return "pair:\n{}\n{}\n".format(self.fw, self.rv)


def load_primers(filename):
    if filename is None or filename == "":
        return []

    f = open(filename, "r")
    primers = []
    i = 0
    for line in f:
        if line[0] == '>' or len(line) < 1:
            continue
        pair = PrimerPair()
        s = line.rstrip().split()
        if line[0] in 'AUTGCN':
            pair.fw.seq = s[0]
            pair.rv.seq = s[1]
            primers.append(pair)
        else:
            primers[-1].fw.left = int(s[0])
            primers[-1].fw.right = int(s[1])
            primers[-1].rv.left = int(s[2])
            primers[-1].rv.right = int(s[3])
            primers[-1].index = i
            i += 1
    return primers


class DataColumn:
    def __init__(self, index):
        self.index = index
        self.array = []
    def __str__(self):
        return "index {}\narray: {}".format(self.index, self.array)


def load_table(filename):
    f = open(filename, 'r')
    headers = f.readline().rstrip().split('\t')
    seq = []
    seq_index = headers.index("Sequence")
    depths = {sample:{cat:OrderedDict() for cat in categories} for sample in samples}
    n_primer_pairs = 0
    seq_len = 0

    primer_reg = re.compile('(Untreated|Modified|Denatured)_primer_pair_[0-9]+_mapped_depth')
    simple_reg = re.compile('(Untreated|Modified|Denatured)_mapped_depth')

    for i, h in enumerate(headers):
        s = h.split('_')
        if h.endswith("_low_mapq_mapped_depth"):
            sample = s[0]
            depths[sample]["excluded"]["low MAPQ"] = DataColumn(i)
        elif h.endswith("_off_target_mapped_depth"):
            sample = s[0]
            depths[sample]["excluded"]["off-target"] = DataColumn(i)
        else:
            m = primer_reg.fullmatch(h)
            if m is not None:
                sample = s[0]
                primer_index = int(s[3])-1
                name = "Primer pair {}".format(primer_index+1)
                depths[sample]["included"][name] = DataColumn(i)
            m = simple_reg.fullmatch(h)
            if m is not None:
                sample = s[0]
                name = "Included"
                depths[sample]["included"][name] = DataColumn(i)

    if False:
        from pprint import pprint
        pprint(depths)
        exit(1)
    
    for line in f:
        s = line.rstrip('\n\r').split('\t')
        seq.append(s[seq_index])
        seq_len += 1
        for cat in categories:
            for sample in samples:
                for name in depths[sample][cat]:
                    i = depths[sample][cat][name].index
                    depths[sample][cat][name].array.append(int(s[i]))

    for cat in categories:
        for sample in samples:
            for name in depths[sample][cat]:
                depths[sample][cat][name].array = np.array(depths[sample][cat][name].array)

    if False:
        from pprint import pprint
        pprint(depths)
        exit(1)

    present_samples = []
    for sample in samples:
        m = 0
        for cat in categories:
            for name in depths[sample][cat]:
                m += np.sum(depths[sample][cat][name].array)
        if m > 0:
            present_samples.append(sample)
    return seq, depths, n_primer_pairs, seq_len, present_samples


def plot_primers(primer_pairs, ax=None):
    if primer_pairs is None:
        return
    if ax is None:
        ax = plt.gca()

    interval = 0.02
    c = 0
    trans = mp.transforms.blended_transform_factory(
                ax.transData, ax.transAxes)
    unsorted = primer_pairs
    sorted_primer_pairs = []

    while len(unsorted) > 0:
        left = 99999999
        ind = None
        for i, p in enumerate(unsorted):
            if p.fw.left < left:
                ind = i
                left = p.fw.left
        sorted_primer_pairs.append(unsorted[ind])
        unsorted = unsorted[:ind] + unsorted[ind+1:]


    for i, p in enumerate(sorted_primer_pairs):
        y = -0.175 - interval * c
        c += 1
        if c % 3 == 0:
            c = 0
        ax.arrow(p.fw.left, y,
                 p.rv.right - p.fw.left + 1,
                 0,
                 clip_on=False,
                 color="gray",
                 transform=trans,
                 head_width=0, head_length=0)
        ax.arrow(p.fw.left, y,
                 p.fw.right-p.fw.left+1,
                 0,
                 clip_on=False,
                 color="darkred",
                 transform=trans,
                 head_width=0, head_length=0,
                 lw=1.5)
        ax.arrow(p.rv.right, y,
                 -(p.rv.right-p.rv.left+1),
                 0,
                 clip_on=False,
                 color="darkblue",
                 transform=trans,
                 head_width=0, head_length=0,
                 lw=1.5)
        if len(sorted_primer_pairs) > 1:
            t = ax.text(p.fw.left+(p.rv.right-p.fw.left)/2.0,
                    y,
                    str(p.index+1),
                    fontsize=4, ha="center", va="center",
                    transform=trans,
                    bbox=dict(facecolor='white', lw=0, pad=0.1))




def plot_figure(seq, depths, n_primer_pairs, seq_len, rna_name, primer_pairs, present_samples):
    #fig = plt.figure(figsize=(10,8))
    #axes = fig.subplots(nrows=3,ncols=1, sharey=True)
    fig, axes = plt.subplots(3, 1, sharey=True, figsize=(10,8))
    axes = axes.flatten()
    #plt.suptitle("{} mapped depths\n(does not include post-alignment basecall quality filter)".format(rna_name),
    #             ha="left")
    plt.subplots_adjust(hspace=0.4, left=0.1, right=0.86, bottom=0.1)

    s = "{} mapped depths\n(does not include post-alignment basecall quality filter)".format(rna_name)
    axes[0].text(0.02, 0.95, s,
             fontsize=14,
             transform=fig.transFigure)

    for i, sample in enumerate(samples):
        if sample not in present_samples:
            continue
        ax = axes[i]
        ax.set_title(sample, ha="left", loc="left")
        xlim = (1, seq_len+1)
        ax.set_xlim(xlim)
        for cat in categories:
            for n, name in enumerate(depths[sample][cat]):
                kwargs = {}
                if name == 'off-target':
                    kwargs = dict(color='black', ls='--')
                elif name == "low MAPQ":
                    kwargs = dict(color='red', ls=':')
                d = depths[sample][cat][name].array
                p = ax.plot(list( range(1, seq_len+1) ),
                           d,
                           label=name,
                           **kwargs)
                if cat == "included" and name != "Included":
                    # shade right and left of primer region
                    if primer_pairs is not None:
                        color = p[0].get_color()
                        p = primer_pairs[n]
                        x = np.array(list( range(1, seq_len+1) ))
                        '''
                        ax.fill_between(x,
                                        0,
                                        d,
                                        where=(x<=p.fw.right),
                                        color=color, alpha=0.2)
                        ax.fill_between(x,
                                        0,
                                        d,
                                        where=(x>=p.rv.left),
                                        color=color, alpha=0.2)
                        '''
                        where = (p.fw.right < x) * (x < p.rv.left)
                        if False:
                            if sample == "Modified":
                                print(name)
                                print(np.where(where)[0])
                                print(x[where])
                                print(d[where])
                        ax.fill_between(x[where],
                                    0,
                                    d[where],
                                    color=color, alpha=0.2, lw=0, zorder=1)

    axes[0].legend(loc=9, ncol=1,
                   bbox_to_anchor=(0.93, 0.91),
                   bbox_transform=fig.transFigure)
    for i, sample in enumerate(samples):
        ax = axes[i]
        if sample not in present_samples:
            ax.remove()
            continue

        ymin, ymax = ax.get_ylim()
        ax.set_ylim((0, ymax))
        ax.set_xlim((1, seq_len+1))
        #ax.set_xlabel("Nucleotide")
        ax.set_ylabel("Read depth")
        ax.yaxis.grid(True)
        ax.set_axisbelow(True)
        ax.xaxis.set_ticks_position('bottom') # no xtick.top rcParam in mpl<2.0
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_minor_locator(mp.ticker.AutoMinorLocator())
        yticks = [int(y) for y in ax.get_yticks()]
        formatted_ticks = []
        for val in yticks:
            formatted_ticks.append(metric_abbreviate(val))
        ax.set_yticklabels(formatted_ticks)
        plot_primers(primer_pairs, ax=ax)

    for i, ax in enumerate(axes[::-1]):
        if samples[::-1][i] in present_samples:
            ax.set_xlabel("Nucleotide", labelpad=10)
            break

def format_abundances(depths, primer_pairs, present_samples, fmt="tab"):
    s = ""

    if fmt == "tab":
        s += "Maximum read depth "
    elif fmt == "human":
        s += "Approximate maximum read depth "

    s += "for reads mapping near each amplicon primer pair\n"
    s += "(See *_mapped_depths.pdf for plots)\n"
    if fmt == "human":
        s += "(See *_per-amplicon_abundance.txt for exact numbers)\n"

    if fmt == "tab":
        col_names = [
            "Name",
            "fw_primer",
            "fw_primer_left",
            "fw_primer_right",
            "rv_primer",
            "rv_primer_left",
            "rv_primer_right",
            "max_depth",
        ]
    elif fmt == "human":
        col_names = [
            "pair",
            "forward",
            "left",
            "right",
            "reverse",
            "left",
            "right",
            "depth",
        ]

    col_widths = [
        5,
        30,
        8,
        8,
        30,
        8,
        8,
        7,
    ]

    for sample in samples:
        if sample not in present_samples:
            continue
        s += "\n"+sample+":\n"
        if fmt == "tab":
            header = '\t'.join(col_names)+"\n"
        elif fmt == "human":
            header = '  '+''.join(['{{:<{}}}'.format(w) for w in col_widths]).format(*col_names)
            fields = ['-'*(w-1) for w in col_widths]
            underline = '  '+''.join(['{{:<{}}}'.format(w) for w in col_widths]).format(
                *fields
            )
            header += "\n"+underline+"\n"

        s += header
        for i, name in enumerate(depths[sample]["included"]):
            p = primer_pairs[i]
            d = depths[sample]["included"][name].array
            maxd = np.max(d)
            name = name.replace(' ', '_')
            lf = '\t'.join(['{}']*8) + '\n'
            vals = [
                name,
                p.fw.seq,
                p.fw.left,
                p.fw.right,
                p.rv.seq,
                p.rv.left,
                p.rv.right,
                maxd
            ]
            if fmt == "human":
                vals[-1] = metric_abbreviate(vals[-1])
                vals[0] = vals[0].lstrip("Primer_pair_")
            if fmt == "tab":
                s += lf.format(*vals)
            elif fmt == "human":
                s += '  '+''.join(['{{:<{}}}'.format(w) for w in col_widths]).format(*vals)+"\n"

    return s

if __name__=="__main__":
    if False:
        print(sys.version)
        import matplotlib as mp
        print(mp.__version__)
        exit(1)


    ap = AP()
    ap.add_argument("--rna-name", type=str, required=True)
    ap.add_argument("--tsv", type=str, required=True,
                    help="Tab-delimited input file produced by make_reactivity_profiles.py (format described in docs/file_formats.md")
    ap.add_argument("--primer-locations", type=str, required=False, default="")
    ap.add_argument("--out", type=str, required=True, help="filename for output figure")
    ap.add_argument("--estimated-abundances", type=str, required=False,
                    help="filename for estimated per-primer-pair abundance output table",
                    default="")

    pa = ap.parse_args(sys.argv[1:])

    title = "{} mapped read depths \n(Does not include post-alignment basecall quality filter)"
    title = title.format(pa.rna_name)

    primer_pairs = None
    if pa.primer_locations != "":
        primer_pairs = load_primers(pa.primer_locations)

    seq, depths, n_primer_pairs, seq_len, present_samples = load_table(pa.tsv)

    if pa.estimated_abundances != "":
        open(pa.estimated_abundances, "w").write(format_abundances(depths, primer_pairs, present_samples))
        print(format_abundances(depths, primer_pairs, present_samples, fmt="human"))

    plot_figure(seq, depths, n_primer_pairs, seq_len, pa.rna_name, primer_pairs, present_samples)
    plt.savefig(pa.out)
