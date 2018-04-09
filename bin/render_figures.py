"""
Given:
 summary tab delimited file

Render:
 - simple SHAPE profile
 - shape profile with error bars
 - modified, untreated, denatured mutation rate profiles and stderr shading
 - depths for each sample FIXME: get log-scaled depths working again
 - (above profiles should only be rendered for RNAs
    below some reasonable length to avoid out-of-memory errors.
    In the future, render zoomed out median profiles)
 - depth histograms
 - reactivity, stderr histograms

Also do quality control checks:
 - good depths over some fraction of nucs
 - good mutation rates above background
 - not too many high mutation rates in untreated sample
 - Write any warnings to stdout, and display prominently in figures.
 - In the future, also check for PCR artifacts (clusters of high
   background positions explained by partial local self-complementarity).

"""

# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys, os, argparse
from numpy import isnan, nan, sqrt
from numpy import nanpercentile as percentile
import numpy as np

from math import ceil

import matplotlib as mp
mp.rcParams["font.sans-serif"].insert(0,"Arial")
mp.rcParams["font.family"] = "sans-serif"
mp.rcParams["pdf.fonttype"] = 42 # use TrueType fonts when exporting PDFs
                                 # (embeds most fonts - this is especially
                                 #  useful when opening in Adobe Illustrator)
mp.rcParams['xtick.direction'] = 'out'
mp.rcParams['ytick.direction'] = 'out'
mp.rcParams['legend.fontsize'] = 14
mp.rcParams['grid.color'] = ".8"
mp.rcParams['grid.linestyle'] = '-'
mp.rcParams['grid.linewidth'] = 1
mp.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

rx_color = "red"
bg_color = "blue"
dc_color = "darkgoldenrod"

def metric_abbreviate(num):
    suffixes = {3:'k',
                6:'M',
                9:"G"}
    s = str(num)
    # replace trailing zeros with metric abbreviation
    zero_count = len(s)-len(s.rstrip('0'))
    suffix = ''
    new_string = str(s)
    for num_zeros in sorted(suffixes.keys()):
        if num_zeros <= zero_count:
            suffix = suffixes[num_zeros]
            new_string = s[:-num_zeros]
    new_string = new_string+suffix
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

def qc_stats(seq,
             rx_rate, bg_rate, dc_rate,
             rx_depth, bg_depth, dc_depth,
             min_depth, max_bg,
             min_depth_pass,
             max_high_bg,
             min_positive,
             high_mut_thresh,
             min_high_mut):
    # FIXME: set up some dummy data for common QC scenarios
    
    print("Quality control checks:\n")

    if rx_rate is None and bg_rate is None and dc_rate is None:
        # no coverage at all
        print("FAIL: no reads mapped to this RNA") 
        return False

    depth_mask_pass = np.full((len(rx_depth),), True, dtype=bool)
    rx_depth_pass = np.full((len(rx_depth),), True, dtype=bool)
    bg_depth_pass = np.full((len(rx_depth),), True, dtype=bool)
    dc_depth_pass = np.full((len(rx_depth),), True, dtype=bool)

    for i in range(len(rx_depth)):
        if rx_rate is not None:
            if rx_depth[i] < min_depth:
                depth_mask_pass[i] = False
                rx_depth_pass[i] = False
        if bg_rate is not None:
            if bg_depth[i] < min_depth:
                depth_mask_pass[i] = False
                bg_depth_pass[i] = False
        if dc_rate is not None:
            if dc_depth[i] < min_depth:
                depth_mask_pass[i] = False
                dc_depth_pass[i] = False

    masked = np.array([c.islower() for c in seq])
    depth_mask_pass = depth_mask_pass * (1-masked)

    
    print("(See README.md for explanation)\n")

    overall_pass = True

    num_depth_pass = sum(depth_mask_pass)
    total_unmasked = np.sum(1-masked)
    depth_pass_frac = num_depth_pass/float(total_unmasked)
    rx_pass_frac = sum(rx_depth_pass)/float(total_unmasked)
    bg_pass_frac = sum(bg_depth_pass)/float(total_unmasked)
    dc_pass_frac = sum(dc_depth_pass)/float(total_unmasked)
    print("Read depth check:")
    msg = "{:.1f}% ({}/{}) nucleotides meet the minimum read depth of {}"
    msg = msg.format(depth_pass_frac*100,
                     num_depth_pass,
                     total_unmasked,
                     min_depth)
    print(msg)
    if depth_pass_frac >= min_depth_pass:
        print("PASS\n")
    else:
        overall_pass = False
        msg = "FAIL: Read depths are too low for accurate reactivity profile creation.\n"
        problem_samples = []
        if rx_pass_frac < min_depth_pass:
            problem_samples.append("modified")
        if bg_rate is not None and bg_pass_frac < min_depth_pass:
            problem_samples.append("untreated")
        if dc_rate is not None and dc_pass_frac < min_depth_pass:
            problem_samples.append("denatured")
        if len(problem_samples) == 1:
            msg += "      The {} sample is particulary low.\n".format(problem_samples[0])
        elif len(problem_samples) > 1:
            msg += "      Problematic samples: "+', '.join(problem_samples)+'.\n'
        msg += "      Check alignment stats to see the amount of target sequence present\n"
        msg += "      in each sample. Better target enrichment or recovery and/or\n"
        msg += "      additional sequencing could help solve this problem.\n"
        # TODO: check alignment stats automatically and warn if poor target recovery
        print(msg)

    msg = "      Possible causes:\n"
    msg += "       - DNA contamination\n"
    msg += "       - poor mixing of chemical reagents and RNA and/or poor\n"
    msg += "         reagent diffusion (if modifying in cells), resulting\n"
    msg += "         in low modification rates\n"
    msg += "       - expired reagents, resulting in low modification rates\n"
    msg += "       - poor reverse transcription conditions, resulting in\n"
    msg += "         low adduct read-through\n"
    msg += "       - extremely highly structured RNA\n"
    causes_msg = msg
    printed_causes = False

    if bg_rate is not None and num_depth_pass>10:
        diff = rx_rate-bg_rate
        num_positive = 0
        for i in range(len(rx_rate)):
            if depth_mask_pass[i] and diff[i] > 0:
                num_positive += 1
        positive_frac = num_positive/float(num_depth_pass)
        print("Mutation rate check:")
        msg = "{:.1f}% ({}/{}) nucleotides have positive mutation rates\n"
        msg += "above background"
        msg = msg.format(positive_frac*100,
                         num_positive,
                         num_depth_pass)
        print(msg)
        if positive_frac >= min_positive:
            msg =  "PASS: There is a clear difference in mutation rates between\n"
            msg += "      modified and untreated samples.\n"
        else:
            overall_pass = False
            msg = "FAIL: Mutation rates do not show a clear difference between\n"
            msg += "      modified and untreated samples.\n"
            msg += causes_msg
            print(msg)
            printed_causes = True

        if False:
            msg = '\nptile\trate diff\n'
            good_diff = np.extract(depth_mask_pass, diff)
            for p in np.linspace(0, 100, 20):
                d = percentile(good_diff, p)
                msg += '{}\t{:.6}\n'.format(p, d)
            msg += '\n'
            print(msg)
        #diff_85th = percentile(good_diff, 85)
        #msg = "85th percentile diff rate: {}\n".format(diff_85th)
        #print(msg)

        high_bg = np.full((len(rx_depth),), False, dtype=bool)
        for i in range(len(high_bg)):
            if depth_mask_pass[i] and bg_rate[i] > max_bg:
                high_bg[i] = True
        num_high_bg = sum(high_bg)
        high_bg_frac = num_high_bg/float(num_depth_pass)

        print("High background check:")
        msg = "{:.1f}% ({}/{}) nucleotides have high background\n"
        msg += "mutation rates."
        msg = msg.format(high_bg_frac*100,
                         num_high_bg,
                         num_depth_pass)
        print(msg)
        if high_bg_frac <= max_high_bg:
            print("PASS: Not too many nucleotides with high background mutation rates.\n")
        else:
            overall_pass = False
            msg = "FAIL: Too many nucleotides with high background mutation rates.\n"
            msg += "      This can be caused by the presence of native modifications\n"
            msg += "      or sequence variants.\n"
            print(msg)

    if bg_rate is None:
        # if just a single sample, can at least check that mutation rates
        # are as high as expected.
        diff = rx_rate

    if num_depth_pass>10:
        num_high_mut = 0
        num_inv_high_mut = 0
        # TODO: could do this faster with numpy array boolean ops
        for i in range(len(diff)):
            if depth_mask_pass[i]:
                if diff[i] > high_mut_thresh:
                    num_high_mut += 1
                if diff[i] < -high_mut_thresh:
                    num_inv_high_mut += 1
        high_mut_frac = num_high_mut/float(num_depth_pass)
        inv_high_mut_frac = num_inv_high_mut/float(num_depth_pass)
        print("Number highly reactive check:")
        msg = "{:.1f}% ({}/{}) nucleotides show high apparent reactivity."
        msg = msg.format(high_mut_frac*100,
                         num_high_mut,
                         num_depth_pass)
        print(msg)

        if high_mut_frac >= min_high_mut:
            if inv_high_mut_frac < min_high_mut:
                msg = "PASS: The expected number of nucleotides or more are highly\n"
                msg += "      reactive."
                print(msg)
            else:
                overall_pass = False
                print("FAIL: modified and untreated samples might be swapped\n")
        else:
            overall_pass = False
            if printed_causes:
                print("FAIL: see possible causes listed above\n")
            else:
                print("FAIL\n"+causes_msg)

    return overall_pass

def render_profiles(num, seq, reactivity, stderr,
                    rx_rates, bg_rates, dc_rates,
                    rx_depth, bg_depth, dc_depth,
                    rx_simple_depth, bg_simple_depth, dc_simple_depth,
                    fileout, name, qc_pass):
    # FIXME: this is fairly ugly code - should at least break each panel into a separate func

    if rx_rates is None and bg_rates is None and dc_rates is None:
        no_mapped = True
    else:
        no_mapped = False

    legend_labels = []
    if rx_rates is not None:
        legend_labels.append("Modified")
        rx_err = sqrt(rx_rates) / sqrt(rx_depth)
    else:
        rx_rates = np.zeros((len(rx_depth),))
        rx_err = np.zeros((len(rx_depth),))
    if bg_rates is not None:
        legend_labels.append("Untreated")
        bg_err = sqrt(bg_rates) / sqrt(bg_depth)
    else:
        bg_rates = np.zeros((len(rx_depth),))
        bg_err = np.zeros((len(rx_depth),))
    if dc_rates is not None:
        legend_labels.append("Denatured")
        dc_err = sqrt(dc_rates) / sqrt(dc_depth)
    else:
        dc_rates = np.zeros((len(rx_depth),))
        dc_err = np.zeros((len(rx_depth),))

    # Add a zeroeth nuc so axis numbering works correctly
    # There's probably a better way to do this
    num = np.append(0, num)
    if reactivity is not None:
        reactivity = np.append(0, reactivity)
        stderr = np.append(0, stderr)
    rx_depth = np.append(0, rx_depth)
    bg_depth = np.append(0, bg_depth)
    dc_depth = np.append(0, dc_depth)
    rx_simple_depth = np.append(0, rx_simple_depth)
    bg_simple_depth = np.append(0, bg_simple_depth)
    dc_simple_depth = np.append(0, dc_simple_depth)
    rx_rates = np.append(0, rx_rates)
    bg_rates = np.append(0, bg_rates)
    dc_rates = np.append(0, dc_rates)
    rx_err = np.append(0, rx_err)
    bg_err = np.append(0, bg_err)
    dc_err = np.append(0, dc_err)

    if reactivity is not None:
        orange_thresh = 0.4
        red_thresh = 0.85

        gray_vals = []
        gray_nums = []
        gray_errs = []
        black_vals = []
        black_nums = []
        black_errs = []
        orange_vals = []
        orange_nums = []
        orange_errs = []
        red_vals = []
        red_nums = []
        red_errs = []
        for i in range(len(reactivity)):
            if isnan(reactivity[i]):
                gray_vals.append(-1)
                gray_nums.append(num[i])
                gray_errs.append(0)
            elif reactivity[i] < orange_thresh:
                black_vals.append(reactivity[i])
                black_nums.append(num[i])
                black_errs.append(stderr[i])
            elif reactivity[i] < red_thresh:
                orange_vals.append(reactivity[i])
                orange_nums.append(num[i])
                orange_errs.append(stderr[i])
            else:
                red_vals.append(reactivity[i])
                red_nums.append(num[i])
                red_errs.append(stderr[i])

    yMin, ymax = (-0.5, 4)
    left_inches = 0.9
    right_inches = 0.4
    sp_width = len(num)*0.032
    fig_width = max(7,sp_width+left_inches+right_inches)
    fig = plt.figure(figsize=(fig_width,8))

    left_percent = left_inches/fig_width
    right_percent = 1-right_inches/fig_width
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(313)
    ax3 = plt.subplot(312)
    plt.subplots_adjust(hspace=0.5, left=left_percent,right=right_percent,top=0.94)

    near_black = (0,0,1/255.0)

    if reactivity is not None:
        if len(gray_nums)>0:
            ax1.bar(gray_nums,gray_vals,
                     align="center",
                     width=1.05, color="0.80", edgecolor="0.80",linewidth=0.0)
        if len(black_nums)>0:
            ax1.bar(black_nums,black_vals,
                     align="center",
                     width=1.05, color="black", edgecolor="black",linewidth=0.0,
                     yerr=black_errs,ecolor=near_black,capsize=1)
        if len(orange_nums)>0:
            ax1.bar(orange_nums,orange_vals,
                     align="center",
                     width=1.05, color="orange",edgecolor="orange",linewidth=0.0,
                     yerr=orange_errs,ecolor=near_black,capsize=1)
        if len(red_nums)>0:
            ax1.bar(red_nums,red_vals,
                     align="center",
                     width=1.05,color="red",edgecolor="red",linewidth=0.0,
                     yerr=red_errs,ecolor=near_black,capsize=1)

    #print("title: "+name)
    ax1title = ax1.set_title(name, horizontalalignment="left", fontsize=16)
    x,y = ax1title.get_position()
    ax1title.set_position((0,y))
    ax1.set_ylim(yMin,ymax)
    ax1.set_xlim(1,len(num))
    #ax1.set_yticks(fontsize=9)

    if not qc_pass:
        msg = "Note: possible data quality issue - see log file"
        return_flag = False
        if no_mapped:
            msg = "ERROR: no reads mapped to this RNA"
            return_flag = True
        txt = plt.text(30,573,
                 msg,
                 ha='left', va='top',
                 fontsize=16, color='red',
                 transform=mp.transforms.IdentityTransform())
        if return_flag:
            plt.savefig(fileout)
            return 

    #tickNums = range(num[0]+10,num[-1]+1,10)
    #tickPos = range(num[0]+9,num[-1],10)
    #ax1.set_xticks(tickPos,tickNums,fontsize=9,rotation=30)
    #ax1.set_xticks(fontsize=9)

    ax1.yaxis.grid(True)
    ax1.set_axisbelow(True)

    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    for loc, spine in ax1.spines.items():
        if loc == 'bottom':
            spine.set_position(('outward', 6))  # move outward (down) 6 pts
            spine.set_smart_bounds(True)
    for loc, spine in ax1.spines.items():
        if loc == 'left':
            spine.set_position(('outward', 6))  # move outward (left) 6 pts
            spine.set_smart_bounds(True)

    # need to add labels after moving spines, otherwise they will disappear
    ax1xlabel = ax1.set_xlabel("Nucleotide", horizontalalignment="left", fontsize=14, labelpad=0)
    x,y = ax1xlabel.get_position()
    ax1xlabel.set_position((0,y))
    ax1ylabel = ax1.set_ylabel("Shape Reactivity", horizontalalignment="left", fontsize=14)
    x,y = ax1ylabel.get_position()
    ax1ylabel.set_position((x,0))

    if reactivity is not None:
        # add a SHAPE colorbar to the vertical axis
        # uses a little transformation magic to place correctly
        inv = ax1.transData.inverted()
        for loc, spine in ax1.spines.items():
            if loc == 'left':
                trans = spine.get_transform()
        pt = trans.transform_point([0,0])
        pt2 = inv.transform_point(pt)
        rectX = pt2[0]
        ptA = (0,0)
        ptB = (6,0)
        ptA2 = inv.transform_point(ptA)
        ptB2 = inv.transform_point(ptB)
        rectW = ptB2[0]-ptA2[0]
        rect = Rectangle((rectX,-0.5), rectW, orange_thresh+0.5, facecolor="black", edgecolor="none")
        ax1.add_patch(rect)
        rect.set_clip_on(False)
        rect = Rectangle((rectX,orange_thresh), rectW, red_thresh-orange_thresh, facecolor="orange", edgecolor="none")
        ax1.add_patch(rect)
        rect.set_clip_on(False)
        rect = Rectangle((rectX,red_thresh), rectW, 4-red_thresh, facecolor="red", edgecolor="none")
        ax1.add_patch(rect)
        rect.set_clip_on(False)

    ax1.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax1.get_yaxis().tick_left()

    ax1.tick_params(axis='y',which='minor',left='off')
    #ax1.tick_params(axis='x',which='minor')

    ax1.minorticks_on()

    yticks = ax1.get_yticks()
    stripped_ticks = [str(val).rstrip('0').rstrip('.') for val in yticks]
    ax1.set_yticklabels(stripped_ticks)

    for line in ax1.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax1.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax1.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)


    # put nuc sequence below axis
    font_prop = mp.font_manager.FontProperties(family = "monospace", style="normal",weight="bold",size="4.5")
    for i in range(seq.shape[0]):
        nuc = seq[i]
        if nuc == "T":
            nuc = "U"
        color_dict = {"A": "#f20000",
                      "U": "#f28f00",
                      "G": "#00509d",
                      "C": "#00c200"}
        if nuc in color_dict:
            col = color_dict[nuc]
        elif nuc.upper() in color_dict:
            col = color_dict[nuc.upper()]
        else:
            col = "black"
        ax1.annotate(nuc, xy=(i+1, -0.67),fontproperties=font_prop,color=col,annotation_clip=False, horizontalalignment="center")

    handles = []
    h, = ax2.plot(num, rx_simple_depth, linewidth = 1.5, color=rx_color, alpha=1.0)
    ax2.plot(num, rx_depth, linewidth = 1.0, color=rx_color, alpha=0.3)
    handles.append(h)
    h, = ax2.plot(num, bg_simple_depth, linewidth = 1.5, color=bg_color, alpha=1.0)
    ax2.plot(num, bg_depth, linewidth = 1.0, color=bg_color, alpha=0.3)
    handles.append(h)
    h, = ax2.plot(num, dc_simple_depth, linewidth = 1.5, color=dc_color, alpha=1.0)
    ax2.plot(num, dc_depth, linewidth = 1.0, color=dc_color, alpha=0.3)
    handles.append(h)
    ax2.set_xlim(1,len(num))
    #ax2.legend(["+Reagent","Background","Denatured"], bbox_to_anchor=(1.1,1.1))
    leg = ax2.legend(handles, legend_labels, loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.75)
    xmin, xmax, ymin, ymax = ax2.axis()
    ax2.set_ylim(0,ymax)
    #ax2.set_yscale('log')
    #ax2.set_yscale('symlog')# useful, but disabled because of a matplotlib/pyparsing bug
    ax2xlabel = ax2.set_xlabel("Nucleotide\n(note: effective read depths shown in lighter colors)", horizontalalignment="left", fontsize=14, labelpad=0)
    x,y = ax2xlabel.get_position()
    ax2xlabel.set_position((0,y))

    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax2.get_yaxis().tick_left()

    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='minor',left='off')
    #ax2.tick_params(axis='x',which='minor')

    #xlabels = ["%.2f"%v for v in xticks]
    #ax3.set_xticks(xticks)
    #ax3.set_xticklabels(xlabels,rotation = -45, horizontalalignment='left')

    yticks = [int(y) for y in ax2.get_yticks()]
    formatted_ticks = []
    for val in yticks:
        formatted_ticks.append(metric_abbreviate(val))
    ax2.set_yticklabels(formatted_ticks)

    for line in ax2.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax2.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax2.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    ax2.yaxis.grid(True)
    ax2.set_axisbelow(True)

    ax2ylabel = ax2.set_ylabel("Read depth", horizontalalignment="left", fontsize=14)
    x, y = ax2ylabel.get_position()
    ax2ylabel.set_position((x, 0))
    # tried to make an offset, smaller font note about effective depths,
    # but couldn't get positioning/transforms to work properly.
    # For now just putting in xaxis label

    # choose a decent range for axis, excluding high-background positions
    good_indices = []
    for i in range(len(bg_rates)):
        if bg_rates[i]<=0.05 or isnan(bg_rates[i]):
            good_indices.append(i)
    temp_rates = [rx_rates[i] for i in good_indices]
    near_top_rate = percentile(temp_rates,98.0)
    maxes = [0.32,0.16,0.08,0.04,0.02,0.01]
    ymax = maxes[0]
    for i in range(len(maxes)):
        if near_top_rate<maxes[i]:
            ymax = maxes[i]

    rx_upper = rx_rates + rx_err
    rx_lower = rx_rates - rx_err
    bg_upper = bg_rates + bg_err
    bg_lower = bg_rates - bg_err
    dc_upper = dc_rates + dc_err
    dc_lower = dc_rates - dc_err

    ax3xlabel = ax3.set_xlabel("Nucleotide", horizontalalignment="left", fontsize=14, labelpad=0)
    x,y = ax3xlabel.get_position()
    ax3xlabel.set_position((0,y))
    ax3ylabel = ax3.set_ylabel("Mutation rate (%)", horizontalalignment="left", fontsize=14)
    x,y = ax3ylabel.get_position()
    ax3ylabel.set_position((x,0))

    ax3.plot(num, rx_rates, zorder=3, color=rx_color, linewidth=1.5)
    ax3.plot(num, bg_rates, zorder=2, color=bg_color, linewidth=1.5)
    ax3.plot(num, dc_rates, zorder=2, color=dc_color, linewidth=1.5)
    ax3.fill_between(num, rx_lower, rx_upper, edgecolor="none", alpha=0.5, facecolor=rx_color)
    ax3.fill_between(num, bg_lower, bg_upper, edgecolor="none", alpha=0.5, facecolor=bg_color)
    ax3.fill_between(num, dc_lower, dc_upper, edgecolor="none", alpha=0.5, facecolor=dc_color)
    ax3.legend(legend_labels, loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.75)
    ax3.set_xlim((1,len(rx_rates)))
    ax3.set_ylim((0,ymax))

    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    ax3.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax3.get_yaxis().tick_left()

    ax3.minorticks_on()
    ax3.tick_params(axis='y',which='minor',left='off')

    ticks = [x*100 for x in ax3.get_yticks()]
    ax3.set_yticklabels([str(val).rstrip('0').rstrip('.') for val in ticks])

    for line in ax3.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax3.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax3.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    ax3.yaxis.grid(True)
    ax3.set_axisbelow(True)

    # TODO: add a tick for the first nuc - can't seem to add one without screwing
    # up all the other ticks

    plt.savefig(fileout)


def draw_median(ax, vals, col, int_only=False):
    xmin, xmax, ymin, ymax = ax.axis()
    med = np.nanmedian(vals)
    ax.axvline(x=med, color=col, zorder=0, alpha=0.7)
    txt = str(med)
    if int_only == True:
        txt = metric_abbreviate(int(med))
    ax.annotate(txt, xy=(med, (ymax - ymin) * 0.6), color=col, rotation=-45, fontweight="bold")


def draw_percentile(ax, vals, percent, col, x1=0, x2=0, y=0, percentage=False, name=""):
    xmin, xmax, ymin, ymax = ax.axis()
    # med = calcQuartile(vals,percentile)
    med = percentile(vals, percent)
    if percentage:
        med = med * 100
    if med == int(med):
        txt = "%d" % med
    else:
        txt = "%0.2f" % med
        if percentage:
            txt = txt.rstrip('0')
    if not percentage:
        txt = commify("%i" % med)
    if percentage:
        txt = txt + "%"
    ax.text(x1, y, name, transform=ax.transAxes, fontsize=10, horizontalalignment="left")
    ax.text(x2, y, txt, transform=ax.transAxes, fontsize=10, horizontalalignment="right")


def write_histograms(shape, stderr,
                     min_depth, max_bg,
                     rx_rate, bg_rate, dc_rate,
                     rx_depth, bg_depth, dc_depth,
                     fileout, name, qc_pass):

    # limit rate histograms to positions included in final reactivity profile
    #if rx_rate is not None:
    #    rx_rate = rx_rate[np.isfinite(shape)]
    #if bg_rate is not None:
    #    bg_rate = bg_rate[np.isfinite(shape)]
    #if dc_rate is not None:
    #    dc_rate = dc_rate[np.isfinite(shape)]

    tlabel_size = 10

    legend_labels = []
    if rx_rate is not None:
        legend_labels.append("Modified")
    if bg_rate is not None:
        legend_labels.append("Untreated")
    if dc_rate is not None:
        legend_labels.append("Denatured")

    fig, (ax1,ax2,ax3) = plt.subplots(nrows=1, ncols=3)
    fig.set_size_inches(10,6)
    fig.subplots_adjust(bottom=0.6, top=0.85, wspace=0.5, left=0.08, right=0.95)

    if len(name)<20:
        title = plt.suptitle(name,fontsize=18,horizontalalignment="left", x=0.02)
    if not qc_pass:
        #bounds = title.get_window_extent(fig.canvas.get_renderer()).get_points()
        #print("title bounds: {}".format(bounds))
        if len(name)<20:
            offset = len(name)*12
        else:
            offset = 0
        msg = "Note: possible data quality issue - see log file"
        return_flag = False        
        if bg_rate is None and rx_rate is None and dc_rate is None:
            return_flag = True
            msg = "ERROR: no reads mapped to this RNA"
        txt = plt.text(15+offset, 430,
                 msg,
                 ha='left', va='top',
                 fontsize=16, color='red',
                 transform=mp.transforms.IdentityTransform())

        #bounds = txt.get_window_extent(fig.canvas.get_renderer()).get_points()
        #print("txt bounds: {}".format(bounds))
        if return_flag:
            plt.savefig(fileout)
            return

    ptiles = []
    if rx_rate is not None and len(rx_rate)>10:
        ptiles.append(percentile(rx_rate, 90.0))
    if bg_rate is not None and len(bg_rate)>10:
        ptiles.append(percentile(bg_rate, 90.0))
    if dc_rate is not None and len(dc_rate)>10:
        ptiles.append(percentile(dc_rate, 90.0))
    max90percentile = max(ptiles)

    num_bins = 30
    # use fewer bins for short RNAs
    if len(rx_depth)<500:
        num_bins = 10
    int_bins = range(0,num_bins+1,1)

    # adjust rate axis for higher mutation rates (e.g. from DMS or other highly mutagenic reagents)
    bin_max = 0.03
    if max90percentile < 0.01:
        bin_max = 0.008
    bin_width = bin_max/(num_bins)
    float_bins = [b*bin_width for b in int_bins]

    ax1.set_ylabel("Normalized\nnucleotide count",fontsize=13)
    ax1.set_xlabel("Mutation rate (%)", fontsize=13)
    #title = ax1.set_title("Mutation rates",fontsize=15)
    ax1.set_title('Mutation rates', x=0.5,y=1.08)
    #y = title.get_y()
    #title.set_y(y

    ax1.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax1.get_yaxis().tick_left()

    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    for line in ax1.get_yticklines() + ax1.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(1)

    ax1.tick_params(axis='both', labelsize=tlabel_size)

    # TODO: simplify these repeated calls

    if "Modified" in legend_labels and len(rx_rate)>10:
        rxpdf, bins, rx_rate_patches = ax1.hist(rx_rate[np.isfinite(rx_rate)], float_bins, histtype='step', lw=1.5,color=rx_color,zorder=4)
        max_rx_rate_y = max([x[1] for x in rx_rate_patches[0].get_xy()])
        # rescale steplots to max out at 1
        rx_rate_patches[0].set_xy([[c[0],c[1]/max_rx_rate_y] for c in rx_rate_patches[0].get_xy()])
    if "Untreated" in legend_labels and len(bg_rate)>10:
        bgpdf,bins, bg_rate_patches = ax1.hist(bg_rate[np.isfinite(bg_rate)], float_bins, histtype='step', lw=1.5, color=bg_color,zorder=3)
        max_bg_rate_y = max([x[1] for x in bg_rate_patches[0].get_xy()])
        bg_rate_patches[0].set_xy([[c[0],c[1]/max_bg_rate_y] for c in bg_rate_patches[0].get_xy()])
    if "Denatured" in legend_labels and len(dc_rate)>10:
        dcpdf,bins, dc_rate_patches = ax1.hist(dc_rate[np.isfinite(dc_rate)], float_bins, histtype='step', lw=1.5, color=dc_color,zorder=2)
        max_dc_rate_y = max([x[1] for x in dc_rate_patches[0].get_xy()])
        dc_rate_patches[0].set_xy([[c[0],c[1]/max_dc_rate_y] for c in dc_rate_patches[0].get_xy()])

    ymax = 1.0
    xmax = bin_max

    leg = ax1.legend(legend_labels, framealpha=0.75, fontsize=11)
    for l in leg.get_lines():
        l.set_linewidth(1.5)
        l.set_alpha(1.0)

    ax1.set_xlim([0.0,xmax])
    ax1.set_ylim([0,ymax])

    ticks = [x*100 for x in ax1.get_xticks()]
    ax1.set_xticklabels(["%d"%val if val==int(val) else "%s"%val for val in ticks])

    if "Modified" in legend_labels and len(rx_rate)>10:
        draw_percentile(ax1, rx_rate, 95.0, "black", x1=0, x2=1, y=-1, percentage=True, name="Modified sample:\n 95th percentile rate:")
        draw_percentile(ax1, rx_rate, 5.0, "black", x1=0, x2=1, y=-1.1, percentage=True, name=" Median rate:")

    if "Untreated" in legend_labels and len(bg_rate)>10:
        draw_percentile(ax1, bg_rate, 95.0, "black", x1=0, x2=1, y=-1.4, percentage=True, name="Untreated sample:\n 95th percentile rate:")
        draw_percentile(ax1, bg_rate, 5.0, "black", x1=0, x2=1, y=-1.5, percentage=True, name=" Median rate:")

    if "Denatured" in legend_labels and len(dc_rate)>10:
        draw_percentile(ax1, dc_rate, 95.0, "black", x1=0, x2=1, y=-1.8, percentage=True, name="Denatured sample:\n 95th percentile rate:")
        draw_percentile(ax1, dc_rate, 5.0, "black", x1=0, x2=1, y=-1.9, percentage=True, name=" Median rate:")

    # TODO: also generate percentile stats on the modified-untreated diff

    ax2.set_title("Read depths", x=0.5,y=1.08)
    ax2.set_ylabel("Normalized\nnucleotide count",fontsize=13)
    ax2.set_xlabel("Effective read depth", fontsize=13)

    max90percentile = max([percentile(bg_depth, 90.0),
                           percentile(rx_depth, 90.0),
                           percentile(dc_depth, 90.0)])

    xmax = max90percentile*1.1

    if "Modified" in legend_labels:
        rxpdf, bins, rx_depth_patches = ax2.hist(rx_depth, bins=num_bins, range=(0, xmax), histtype='step', lw=1.5,
                                                 color=rx_color, zorder=4)
        max_rx_depth_y = max([x[1] for x in rx_depth_patches[0].get_xy()])
        # rescale steplots to max out at 1
        rx_depth_patches[0].set_xy([[c[0],c[1]/max_rx_depth_y] for c in rx_depth_patches[0].get_xy()])
    if "Untreated" in legend_labels:
        bgpdf, bins, bg_depth_patches = ax2.hist(bg_depth, bins=num_bins, range=(0, xmax), histtype='step', lw=1.5,
                                                 color=bg_color, zorder=3)
        max_bg_depth_y = max([x[1] for x in bg_depth_patches[0].get_xy()])
        bg_depth_patches[0].set_xy([[c[0],c[1]/max_bg_depth_y] for c in bg_depth_patches[0].get_xy()])
    if "Denatured" in legend_labels:
        dcpdf, bins, dc_depth_patches = ax2.hist(dc_depth, bins=num_bins, range=(0, xmax), histtype='step', lw=1.5,
                                                 color=dc_color, zorder=2)
        max_dc_depth_y = max([x[1] for x in dc_depth_patches[0].get_xy()])
        dc_depth_patches[0].set_xy([[c[0],c[1]/max_dc_depth_y] for c in dc_depth_patches[0].get_xy()])

    ymax = 1.0

    #leg = ax2.legend(["Modified","Untreated","Denatured"], framealpha=0.75, fontsize=11)
    #for l in leg.get_lines():
    #    l.set_linewidth(1.5)
    #    l.set_alpha(1.0)

    ax2.set_xlim([0.0,xmax])
    ax2.set_ylim([0.0,ymax])

    if "Modified" in legend_labels:
        draw_percentile(ax2, rx_depth, 50.0, "black", x1=0, x2=1, y=-1, name="Modified sample:\n Median depth:")
        draw_percentile(ax2, rx_depth, 5.0, "black", x1=0, x2=1, y=-1.1, name=" 5th percentile depth:")
    if "Untreated" in legend_labels:
        draw_percentile(ax2, bg_depth, 50.0, "black", x1=0, x2=1, y=-1.4, name="Untreated sample:\n Median depth:")
        draw_percentile(ax2, bg_depth, 5.0, "black", x1=0, x2=1, y=-1.5, name=" 5th percentile depth:")
    if "Denatured" in legend_labels:
        draw_percentile(ax2, dc_depth, 50.0, "black", x1=0, x2=1, y=-1.8, name="Denatured sample:\n Median depth:")
        draw_percentile(ax2, dc_depth, 5.0, "black", x1=0, x2=1, y=-1.9, name=" 5th percentile depth:")



    ax2.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax2.get_yaxis().tick_left()

    ax2.tick_params(axis='both', labelsize=tlabel_size)

    xticks = [int(x) for x in ax2.get_xticks()]
    formatted_ticks = []
    for val in xticks:
        formatted_ticks.append(metric_abbreviate(val))
    ax2.set_xticklabels(formatted_ticks, rotation=90)

    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    for line in ax2.get_yticklines() + ax2.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(1)

    ax3.set_xlabel("Normalized reactivity", fontsize=13)
    ax3.set_ylabel("Normalized\nnucleotide count", fontsize=13)
    ax3.set_title("Reactivity distribution", x=0.5,y=1.08)
    #xticks = [x/10.0 for x in range(-5, 30, 5)] + [0.4, 0.85]
    #xticks = [-1.0, -0.5, 0.0, 0.4, 0.85, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    xticks = [-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
    xlabels = ["%d"%val if val==int(val) else "%s"%val for val in xticks]

    ax3.tick_params(axis='both', labelsize=tlabel_size)

    ax3.set_xticks(xticks)
    ax3.set_xticklabels(xlabels)

    if shape is not None:
        pdf, bins,shape_patches = ax3.hist(shape[np.isfinite(shape)], bins=num_bins, range=(-1,4), histtype='step', lw=1.5, color="black", zorder=4)
        pdf, bins, err_patches = ax3.hist(stderr[np.isfinite(shape)], bins=num_bins, range=(-1,4), histtype='step', lw=1.5, color='blue', zorder=3)

        max_shape_y = max([x[1] for x in shape_patches[0].get_xy()])
        max_err_y = max([x[1] for x in err_patches[0].get_xy()])

        # rescale steplots to max out at 1
        shape_patches[0].set_xy([[c[0],c[1]/max_shape_y] for c in shape_patches[0].get_xy()])
        err_patches[0].set_xy([[c[0],c[1]/max_err_y] for c in err_patches[0].get_xy()])

        leg = ax3.legend(["Reactivities","Stderrs"], framealpha=0.75, fontsize=11)

        ax3.axvline(x=0,ls='-',lw=1,color="0.5",zorder=0)
        #ax3.axvline(x=0.4,ls='-',lw=1,color="0.5",zorder=0)
        #ax3.axvline(x=0.85,ls='-',lw=1,color="0.5",zorder=0)

    #ymin, ymax = ax3.get_ylim()
    #ymax = max(shapePdf)
    ymax = 1
    #ax3.set_xlim([-0.5,3])
    ax3.set_ylim([0,ymax])

    ax3.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax3.get_yaxis().tick_left()

    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    for line in ax3.get_yticklines() + ax3.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(1)

    plt.savefig(fileout)

    
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
    expected_headers_types = [
        ('Nucleotide','int'),
        ('Sequence','str'),

        ('Modified_mutations','int'),
        ('Modified_read_depth','int'),
        ('Modified_effective_depth','int'),
        ('Modified_rate','float'),

        ('Untreated_mutations','int'),
        ('Untreated_read_depth','int'),
        ('Untreated_effective_depth','int'),
        ('Untreated_rate','float'),

        ('Denatured_mutations','int'),
        ('Denatured_read_depth','int'),
        ('Denatured_effective_depth','int'),
        ('Denatured_rate','float'),

        ('Reactivity_profile','float'),
        ('Std_err','float'),

        ('HQ_profile','float'),
        ('HQ_stderr','float'),

        ('Norm_profile','float'),
        ('Norm_stderr','float'),
    ]
    expected_headers, types = zip(*expected_headers_types)
    for x in expected_headers[:-2]:
        if x not in headers:
            raise RuntimeError("File \""+filename+"\" does not contain the expected columns.")

    # match expected columns with actual header in case
    # file format has additional columns or rearrangements
    remapped_types = ['str' for x in range(len(headers))]
    for n in range(len(headers)):
        try:
            i = expected_headers.index(headers[n])
            remapped_types[n] = types[i]
        except ValueError:
            pass

    data = {headers[j]:np.empty(length,dtype=remapped_types[j]) for j in range(len(headers))}
    #sys.stderr.write("data['Modified_occluded_depth']: {}\n".format(data["Modified_occluded_depth"]))

    for i in range(length):
        line = f.readline()
        s = line.strip().split('\t')
        for j in range(len(headers)):
            try:
                #sys.stderr.write("Attempting to remap data array for header {}\n".format(headers[j]))
                #sys.stderr.write("remapped_types[j]: {}\n".format(remapped_types[j]))
                data[headers[j]][i] = eval(remapped_types[j])(s[j])
            #except ValueError:
            #    data[headers[j]][i] = -999 # internally means "no data"
            except IndexError:
                s = "File \""+filename+"\" is misformatted "
                s += "(contains fewer columns than headers)."
                raise RuntimeError(s)
            if i == length-1 and remapped_types[j] == "float":
                if np.isnan(data[headers[j]]).all():
                    data[headers[j]] = None

    return data

if __name__=="__main__":
    parser = argparse.ArgumentParser()

    h = "Tab-delimited file summarizing reactivity"
    h += " profile, mutation rates, and sequencing depth."
    parser.add_argument("--infile", help=h, required=True, type=str)

    h = "PDF file to render reactivity profile, mutation rates, and depths."
    h += " This will only be rendered if the sequence is less than maxlen nucleotides long."
    parser.add_argument("--plot", help=h, type=str)

    h = "Maximum sequence length for rendering profile figures to PDF (large sequences "
    h += "can cause out-of-memory errors)."
    parser.add_argument("--maxlen", help=h, type=int, default=10000)

    h = "PDF file to render summary histograms."
    parser.add_argument("--hist", help=h, type=str)

    h = "Title for figures."
    parser.add_argument("--title", help=h, type=str)

    h = "Minimum required sequencing depth for all provided samples for"
    h += " including a nucleotide position."
    parser.add_argument("--mindepth", help=h, type=int, default=5000)

    h = "Maximum allowed mutation rate in untreated sample (if present)"
    h += " for including a nucleotide position."
    parser.add_argument("--maxbg", help=h, type=float, default=0.05)

    h = "Minimum fraction of nucleotides that must pass mindepth threshold."
    parser.add_argument("--min-depth-pass", help=h, type=float, default=0.8)

    h = "Maximum fraction of nucleotides that can have background mutation rates above maxbg."
    parser.add_argument("--max-high-bg", help=h, type=float, default=0.05)

    h = "Minimum fraction of nucleotides with positive mutation rates above background"
    h += " (excluding low-depth positions)."
    parser.add_argument("--min-positive", help=h, type=float, default=0.5)

    h = "Mutation rate (or rate difference, if using untreated control) threshold"
    h += " for calling highly mutated nucleotides."
    parser.add_argument("--high-mut-thresh", help=h, type=float, default=0.006) # 0.006-0.007

    h = "Minimum fraction of nucleotides that are highly mutated"
    h += " (excluding low-depth positions)."
    parser.add_argument("--min-high-mut", help=h, type=float, default=0.08)

    p = parser.parse_args(sys.argv[1:])

    d = load_tab(p.infile)

    # create output directories if needed
    if p.plot is not None:
        o = os.path.split(p.plot)[0]
        if len(o)>0:
            os.makedirs(o, exist_ok=True)
    if p.hist is not None:
        o = os.path.split(p.hist)[0]
        if len(o)>0:
            os.makedirs(o, exist_ok=True)


    if "Norm_profile" in d:
        profile = d["Norm_profile"]
        stderr = d["Norm_stderr"]
    else:
        profile = None
        stderr = None

    if p.title is not None:
        title = p.title
    else:
        title = ""

    qc_pass = qc_stats(d["Sequence"],
                       d["Modified_rate"],
                       d["Untreated_rate"],
                       d["Denatured_rate"],
                       d["Modified_effective_depth"],
                       d["Untreated_effective_depth"],
                       d["Denatured_effective_depth"],
                       p.mindepth, p.maxbg,
                       p.min_depth_pass,
                       p.max_high_bg,
                       p.min_positive,
                       p.high_mut_thresh,
                       p.min_high_mut)

    if p.plot is not None and (d["Modified_rate"] is None or 
                               len(d["Modified_rate"]) <= p.maxlen):
        render_profiles(
            d["Nucleotide"], d["Sequence"], profile, stderr,
            d["Modified_rate"], d["Untreated_rate"], d["Denatured_rate"],
            d["Modified_effective_depth"], d["Untreated_effective_depth"], d["Denatured_effective_depth"],
            d["Modified_read_depth"], d["Untreated_read_depth"], d["Denatured_read_depth"],
            p.plot, title, qc_pass)

    if p.hist is not None:
        write_histograms(
            profile, stderr,
            p.mindepth, p.maxbg,
            d["Modified_rate"], d["Untreated_rate"], d["Denatured_rate"],
            d["Modified_effective_depth"], d["Untreated_effective_depth"], d["Denatured_effective_depth"],
            p.hist, title, qc_pass)
