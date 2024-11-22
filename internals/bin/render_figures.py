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
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys, os, argparse
from numpy import isnan, nan, sqrt
from numpy import nanpercentile as percentile
import numpy as np

from math import ceil

import matplotlib as mp
mp.use('Agg')
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


class Primer:
    def __init__(self):
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
    for line in f:
        if line[0] == '>' or len(line) < 1:
            continue
        pair = PrimerPair()
        s = line.rstrip().split()
        if line[0] in 'AUTGCN':
            pair.fw.seq = s[0]
            pair.rv.seq = s[1][::-1] # store reverse primer sequence reversed to simplify rendering
            primers.append(pair)
        else:
            primers[-1].fw.left = int(s[0])
            primers[-1].fw.right = int(s[1])
            primers[-1].rv.left = int(s[2])
            primers[-1].rv.right = int(s[3])
    return primers


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
    

    print("SHAPE mode quality control checks:\n")

    if rx_rate is None and bg_rate is None and dc_rate is None:
        # no coverage at all
        print("FAIL: no reads mapped to this RNA") 
        return False   
    
    # make sure a numpy
    seq = np.array(list(seq))


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
            print(msg)
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

        high_bg = np.full((len(rx_depth),), False, dtype=bool)
        for i in range(len(high_bg)):
            if depth_mask_pass[i] and bg_rate[i] > max_bg:
                high_bg[i] = True
        num_high_bg = sum(high_bg)
        high_bg_frac = num_high_bg/float(num_depth_pass)

        print("High background check:")
        msg = "{:.1f}% ({}/{}) nucleotides have high background mutation rates."
        print(msg.format(high_bg_frac*100, num_high_bg, num_depth_pass))

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
            
            # check that reactivities are roughly similar for all nts; if not warn about DMS
            mask = (depth_mask_pass==1) & ( (seq == 'A') | (seq == 'C') ) & np.isfinite(diff)
            ac = np.percentile(diff[mask], 95)
            
            mask = (depth_mask_pass==1) & ( (seq == 'G') | (seq == 'U') ) & np.isfinite(diff)
            gu = np.percentile(diff[mask], 95)

            if ac/gu > 2:
                overall_pass = False
                msg = "WARNING: A/C nts are much more reactive than G/U\n"
                msg += "         This may be a DMS experiment, in which case\n"
                msg += "         DMS mode is recommended"
                print(msg)

            elif inv_high_mut_frac < min_high_mut:
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




def qc_stats_dms(seq,
                rx_rate, bg_rate, dc_rate,
                rx_depth, bg_depth, dc_depth,
                n_rate,
                min_depth, max_bg,
                min_depth_pass,
                max_high_bg,
                min_positive,
                min_high_mut):
    
    high_mut_thresh = {'A':0.05, 'C':0.05, 'U':0.01, 'G':0.002}


    print("\nDMS mode quality control checks:")
    print("--------------------------------")

    if rx_rate is None and bg_rate is not None:
        # Only bg coverage
        print("FAIL: no reads from modified treatment mapped to this RNA")
        return False
    if rx_rate is None and bg_rate is None and dc_rate is None:
        # no coverage at all
        print("FAIL: no reads mapped to this RNA") 
        return False

    # make sure a nump
    seq = np.array(list(seq))

    if bg_rate is not None:
        diff = rx_rate-bg_rate
    else:
        diff = rx_rate
    
    overall_pass = {}
    
    for nt in ('A','C','U','G'):
     
        print("\n***************")
        print("*   Nt = {}    *".format(nt))
        print("***************\n")

        overall_pass[nt] = True

        rx_depth_pass = np.full((len(rx_depth),), True, dtype=bool)
        bg_depth_pass = np.full((len(rx_depth),), True, dtype=bool)
        dc_depth_pass = np.full((len(rx_depth),), True, dtype=bool)

        if rx_rate is not None:
            rx_depth_pass = (rx_depth >= min_depth)
        if bg_rate is not None:
            bg_depth_pass = (bg_depth >= min_depth)
        if dc_rate is not None:
            dc_depth_pass = (dc_depth >= min_depth)
        
        overall_mask = rx_depth_pass & bg_depth_pass & dc_depth_pass & np.char.isupper(seq) & (seq==nt)
        
        
        # check that sufficient number of positions are definied

        num_depth_pass = sum(overall_mask)
        total_unmasked = np.sum(np.char.isupper(seq) & (seq==nt))
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
            overall_pass[nt] = False
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
            print(msg)


        # check high background positions
        if bg_rate is not None and num_depth_pass > 20:
        
            num_high_bg = sum(bg_rate[overall_mask] > max_bg)
            high_bg_frac = num_high_bg/float(num_depth_pass)

            print("High background check:")
            msg = "{:.1f}% ({}/{}) nucleotides have high background mutation rates (>{}) "
            msg += "mutation rates."
            msg = msg.format(high_bg_frac*100,
                             num_high_bg,
                             num_depth_pass, max_bg)
            print(msg)
            if high_bg_frac <= max_high_bg:
                print("PASS: Not too many nucleotides with high background mutation rates.\n")
            else:
                overall_pass[nt] = False
                msg = "FAIL: Too many nucleotides with high background mutation rates.\n"
                msg += "      This can be caused by the presence of native modifications\n"
                msg += "      or sequence variants.\n"
                print(msg)
    

            num_positive = np.sum(diff[overall_mask]>0.0001)
            positive_frac = num_positive/float(num_depth_pass)
            print("Mutation rate check:")
            msg = "{:.1f}% ({}/{}) nucleotides have mutation rates above background"
            msg = msg.format(positive_frac*100,
                             num_positive,
                             num_depth_pass)
            print(msg)
            if positive_frac >= min_positive:
                msg =  "PASS: There is a clear difference in mutation rates between\n"
                msg += "      modified and untreated samples.\n"
                print(msg)
            else:
                overall_pass[nt] = False
                msg = "FAIL: Mutation rates do not show a clear difference between\n"
                msg += "      modified and untreated samples.\n"
                print(msg)

        if num_depth_pass > 20:
            
            num_high_mut = np.sum(diff[overall_mask] > high_mut_thresh[nt])
            num_inv_high_mut = np.sum(diff[overall_mask] < -high_mut_thresh[nt])
            
            high_mut_frac = num_high_mut/float(num_depth_pass)
            inv_high_mut_frac = num_inv_high_mut/float(num_depth_pass)
            
            print("Number highly reactive check:")
            msg = "{:.1f}% ({}/{}) nucleotides have high modification rates (>{})"
            print(msg.format(high_mut_frac*100, num_high_mut, num_depth_pass, high_mut_thresh[nt]))

            if high_mut_frac >= min_high_mut:
                if inv_high_mut_frac < min_high_mut:
                    print("PASS: The expected number of nucleotides (or more) are highly reactive\n")
                else:
                    overall_pass[nt] = False
                    print("FAIL: modified and untreated samples might be swapped\n")
            else:
                overall_pass[nt] = False
                print("FAIL:  Lower than expected reactivity\n")

    
        print("Normalized rate check:")
        mask = np.char.isupper(seq) & (seq == nt)
        if n_rate is None or sum(np.isfinite(n_rate[mask])) == 0:
            msg = "FAIL: No finite normalized rates present.\n"
            msg +=  "      This is either due to the nucleotide failing to\n"
            msg +=  "      pass the normalization factor cutoff of at least\n"
            msg +=  "      0.002 or due to the nucleotide having less than\n"
            msg +=  "      20 nucleotides with quality reactivity information.\n"
            msg +=  "      Please refer to the NormProfile section of this log\n"
            msg +=  "      file for clarification\n"
            overall_pass[nt] = False

        else:
            msg = "PASS: Finite normalized rates present.\n"

        print(msg)

    print("\n-----------------------------------------------------------------")

    print("\n      Nucleotide passes all quality control checks: ")
    for nt in ('A','C','U','G'):
        print("        {}: {}".format(nt, overall_pass[nt]))
    print("\n")

    if (not overall_pass['U'] or not overall_pass['G']) and (overall_pass['A'] or overall_pass['C']):
        msg = "      G and/or U FAILURE, but other reactivities may be usable\n"
        msg += "      Possible causes:\n"
        msg += "       - improper pH control during DMS reaction\n"
        msg += "       - use of a different RT enzyme than MarathonRT\n"
        msg += "       - poor reverse transcription conditions, resulting\n"
        msg += "         in low adduct adduct read-through or detection\n"
        msg += "       - very highly structured or protein bound RNA\n"
        print(msg)


    elif not (overall_pass['A'] and overall_pass['C']):
        msg = "      A and/or C FAILURE: Possible causes:\n"
        msg += "       - low sequencing depth\n"
        msg += "       - DNA contamination\n"
        msg += "       - poor mixing of chemical reagents and RNA and/or poor\n"
        msg += "         reagent diffusion (if modifying in cells), resulting\n"
        msg += "         in low modification rates\n"
        msg += "       - expired reagents, resulting in low modification rates\n"
        msg += "       - poor reverse transcription conditions, resulting in\n"
        msg += "         low adduct read-through\n"
        msg += "       - extremely highly structured RNA\n"
        print(msg)


    return all(overall_pass.values())


def reactivity_graph(ax, name, yMin, ymax, num, no_mapped, reactivity, orange_thresh, red_thresh, seq, qc_pass = None, qc_pass_n7=None, dms = False, N7 = False, fileout=None, message = None, N7_reformat = None):
    if not N7:
       axtitle = ax.set_title(name, horizontalalignment="left", fontsize=16)
       x,y = axtitle.get_position()
       axtitle.set_position((0,y))
    ax.set_ylim(yMin,ymax)
    ax.set_xlim(1,len(num))


    ax.yaxis.grid(True)
    ax.set_axisbelow(True)

    if dms:
      if not N7:
         yticks = np.linspace(0, 1.4, 8)
         yticks = [round(tick, 1) for tick in yticks]
         ax.set_yticks(yticks)
      else:
         yticks = [0, 1.6, 2.3, 3.3]
         ax.set_yticks(yticks)
    else:
       yticks = np.linspace(-.5, 4, 10)
       ax.set_yticks(yticks)


    if not qc_pass and qc_pass != None:
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


    if qc_pass_n7 != None:
       if qc_pass_n7 == "empty":
         print("Normalized N7 data appears to be empty - see log file and ga-profile")
         if N7:
            msg = "Note: Normalized N7 data appears to be empty - see log file and ga-profile"
            txt_n7 = plt.text(60, 350,
                     msg,
                     ha='left', va='top',
                     fontsize=11, color='red',
                     transform=mp.transforms.IdentityTransform())


    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    for loc, spine in ax.spines.items():
        if loc == 'bottom':
            spine.set_position(('outward', 6))  # move outward (down) 6 pts
    for loc, spine in ax.spines.items():
        if loc == 'left':
            spine.set_position(('outward', 6))  # move outward (left) 6 pts

    if dms:
        if not N7:
           if N7_reformat:
              axylabel = ax.set_ylabel("DMS Rctvty(N1,3)", horizontalalignment="left", fontsize=10)
           else:
              axylabel = ax.set_ylabel("DMS Reactivity(N1,3)", horizontalalignment="left", fontsize=14)
        else:
           axylabel = ax.set_ylabel("DMS Rcvty(N7)", horizontalalignment="left", fontsize=10)
    else:
        axylabel = ax.set_ylabel("Shape Reactivity", horizontalalignment="left", fontsize=14)
    x,y = axylabel.get_position()
    axylabel.set_position((x,0))

    if reactivity is not None:
        # add a SHAPE colorbar to the vertical axis
        # uses a little transformation magic to place correctly
        inv = ax.transData.inverted()
        for loc, spine in ax.spines.items():
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
        
        if dms:
            if not N7:
               rect1 = Rectangle((rectX,-0.1), rectW, orange_thresh-yMin, facecolor="black", edgecolor="none")
               rect2 = Rectangle((rectX,orange_thresh), rectW, red_thresh-orange_thresh, facecolor="orange", edgecolor="none")
               rect3 = Rectangle((rectX,red_thresh), rectW, 1.4-red_thresh, facecolor="red", edgecolor="none")

            else:
               rect1 = Rectangle((rectX,-0.22), rectW, orange_thresh-yMin, facecolor="black", edgecolor="none")
               rect2 = Rectangle((rectX,orange_thresh), rectW, red_thresh-orange_thresh, facecolor="pink", edgecolor="none")
               rect3 = Rectangle((rectX,red_thresh), rectW, 3.3-red_thresh, facecolor="purple", edgecolor="none")



        else:
            rect1 = Rectangle((rectX,-0.5), rectW, orange_thresh+0.5, facecolor="black", edgecolor="none")
            rect2 = Rectangle((rectX,orange_thresh), rectW, red_thresh-orange_thresh, facecolor="orange", edgecolor="none")
            rect3 = Rectangle((rectX,red_thresh), rectW, 4-red_thresh, facecolor="red", edgecolor="none")
        
        ax.add_patch(rect1)
        rect1.set_clip_on(False)
        ax.add_patch(rect2)
        rect2.set_clip_on(False)
        ax.add_patch(rect3)
        rect3.set_clip_on(False)


    ax.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax.get_yaxis().tick_left()

    ax.tick_params(axis='y',which='minor',left='off')
    if N7_reformat:
      ax.tick_params(axis='y',which='major', labelsize = 9)

    yticks = ax.get_yticks()
    stripped_ticks = [str(val).rstrip('0').rstrip('.') for val in yticks]
    ax.set_yticklabels(stripped_ticks)

    for line in ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

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

        if dms:
            if not N7:
               ax.annotate(nuc, xy=(i+1, yMin-0.08),fontproperties=font_prop,color=col,annotation_clip=False, horizontalalignment="center")
            else:
               ax.annotate(nuc, xy=(i+1, yMin-0.19),fontproperties=font_prop,color=col,annotation_clip=False, horizontalalignment="center")
        else:
            ax.annotate(nuc, xy=(i+1, -0.67),fontproperties=font_prop,color=col,annotation_clip=False, horizontalalignment="center")


    if message != None:
       warnings = message.split(",")
       if  "low N7" in warnings:
         print("Warning: Unusually low raw N7 reactivity. Consequently, N7 reactivity profiles will not be plotted.\nPlease double check experimental procedure.\n   - This may be caused by highly structured RNA. In this case passing --ignore_low_n7\n    will ignore this QC and plot the N7 data / produce the .mutga file(s) regardless.")
         msg = "Warning: Unusually low raw N7 reactivity. No subplot generated. Please refer to log file."
         txt_n7 = plt.text(60, 347,
                  msg,
                  ha='left', va='top',
                  fontsize=11, color='red',
                  transform=mp.transforms.IdentityTransform())


def mut_graph(ax, dms, num, rx_rates, bg_rates, dc_rates, rx_lower, rx_upper, bg_lower, bg_upper, dc_lower, dc_upper, legend_labels, ymax, N7 = False, N7_reformat = None, message = None):

    mute = False
    if message != None:
       warnings = message.split(",")
       if  "low N7" in warnings:
         mute = True


    if dms:
        if not N7:
           if N7_reformat:
              axylabel = ax.set_ylabel("Mut. Rt.(N1,3)(%)", horizontalalignment="left", fontsize=10)
           else:
              axylabel = ax.set_ylabel("Mutation Rate(N1,3)(%)", horizontalalignment="left", fontsize=14)
        else:
           axylabel = ax.set_ylabel("Mut. Rt.(N7)(%)", horizontalalignment="left", fontsize=10)
    else:
        axylabel = ax.set_ylabel("Mutation rate(%)", horizontalalignment="left", fontsize=14)
    x,y = axylabel.get_position()
    axylabel.set_position((x,0))


    handles = []
    if not N7:
       h, = ax.plot(num, rx_rates, zorder=3, color=rx_color, linewidth=1.5)
       handles.append(h)
       h, = ax.plot(num, bg_rates, zorder=2, color=bg_color, linewidth=1.5)
       handles.append(h)
       h, = ax.plot(num, dc_rates, zorder=2, color=dc_color, linewidth=1.5)
       handles.append(h)

       ax.fill_between(num, rx_lower, rx_upper, edgecolor="none", alpha=0.5, facecolor=rx_color)
       ax.fill_between(num, bg_lower, bg_upper, edgecolor="none", alpha=0.5, facecolor=bg_color)
       ax.fill_between(num, dc_lower, dc_upper, edgecolor="none", alpha=0.5, facecolor=dc_color)
    elif not mute:
       h, = ax.plot(num, rx_rates, zorder=3, color="purple", linewidth=1.5)
       handles.append(h)
       h, = ax.plot(num, bg_rates, zorder=2, color="pink", linewidth=1.5)
       handles.append(h)
       h, = ax.plot(num, dc_rates, zorder=2, color=dc_color, linewidth=1.5)
       handles.append(h)
       ax.fill_between(num, rx_lower, rx_upper, edgecolor="none", alpha=0.5, facecolor="purple")
       ax.fill_between(num, bg_lower, bg_upper, edgecolor="none", alpha=0.5, facecolor="pink")
       ax.fill_between(num, dc_lower, dc_upper, edgecolor="none", alpha=0.5, facecolor=dc_color)
       ax.legend(handles, legend_labels, loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.60)

    ax.set_xlim((1,len(rx_rates)))
    ax.set_ylim((0,ymax))
    if dms and False:
        ax.set_ylim((1e-4, max(max(rx_rates), max(bg_rates))))
        ax.set_yscale('log')

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax.get_yaxis().tick_left()

    ax.tick_params(axis='y',which='minor',left='off')
    if N7_reformat:
       ax.tick_params(axis='y',which='major', labelsize = 9)

    if dms:
      if not N7:
         ax.set_ylim(0, .16)
         ax.set_yticks(np.linspace(0, .16, num = 9))
      else:
         ax.set_ylim(0, .08)
         ax.set_yticks(np.linspace(0, .08, num = 9))

    ticks = [x*100 for x in ax.get_yticks()]
    ticks = [int(tick) for tick in ticks]
    ax.set_yticklabels([str(val) for val in ticks])


    for line in ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    ax.yaxis.grid(True)
    ax.set_axisbelow(True)



def render_profiles(num, seq, reactivity, stderr,
                    rx_rates, bg_rates, dc_rates,
                    rx_depth, bg_depth, dc_depth,
                    rx_simple_depth, bg_simple_depth, dc_simple_depth,
                    fileout, name, qc_pass,
 
                    numN = None , seqN = None , reactivityN = None , stderrN = None ,
                    rx_ratesN = None , bg_ratesN = None , dc_ratesN = None ,
                    rx_depthN = None , bg_depthN = None , dc_depthN = None ,
                    rx_simple_depthN = None , bg_simple_depthN = None , dc_simple_depthN = None ,

                    dms=False, N7=None, qc_pass_n7=None, message = None):
    # FIXME: this is fairly ugly code - should at least break each panel into a separate func
    # TODO: Route all the N stuff through all the pre - processing


    if N7 != None:
      N7_reformat = True

    if rx_rates is None and bg_rates is None and dc_rates is None:
        no_mapped = True
    else:
        no_mapped = False

    legend_labels = []
    if rx_rates is not None:
        legend_labels.append("Modified")
        rx_err = sqrt(rx_rates) / sqrt(rx_depth)
        if N7 != None:
           rx_errN = sqrt(rx_ratesN) / sqrt(rx_depthN)
    else:
        rx_rates = np.zeros((len(rx_depth),))
        if N7 != None:
           rx_ratesN = np.zeros((len(rx_depthN),))
        rx_err = np.zeros((len(rx_depth),))
        if N7 != None:
           rx_errN = np.zeros((len(rx_depthN),))
    if bg_rates is not None:
        legend_labels.append("Untreated")
        bg_err = sqrt(bg_rates) / sqrt(bg_depth)
        if N7 != None:
           bg_errN = sqrt(bg_ratesN) / sqrt(bg_depthN)
    else:
        bg_rates = np.zeros((len(rx_depth),))
        if N7 != None:
           bg_ratesN = np.zeros((len(rx_depthN),))
        bg_err = np.zeros((len(rx_depth),))
        if N7 != None:
           bg_errN = np.zeros((len(rx_depthN),))
    if dc_rates is not None:
        legend_labels.append("Denatured")
        dc_err = sqrt(dc_rates) / sqrt(dc_depth)
        if N7 != None:
           dc_errN = sqrt(dc_ratesN) / sqrt(dc_depthN)
    else:
        dc_rates = np.zeros((len(rx_depth),))
        if N7 != None:
           dc_ratesN = np.zeros((len(rx_depthN),))
        dc_err = np.zeros((len(rx_depth),))
        if N7 != None:
           dc_errN = np.zeros((len(rx_depthN),))

    # Add a zeroeth nuc so axis numbering works correctly
    # There's probably a better way to do this
    num = np.append(0, num)
    #print("This is the zeroeth nucleotide addition step. Num: {}".format(num))
    if N7 != None:
       numN = np.append(0, numN)

    if reactivity is not None:
        reactivity = np.append(0, reactivity)
        stderr = np.append(0, stderr)
        if N7 != None:
           reactivityN = np.append(0, reactivityN)
           stderrN = np.append(0, stderrN)
    if reactivity is None:
        print("Reactivity is none.")



    #TODO: Combine this all into one conditional...
    rx_depth = np.append(0, rx_depth)
    if N7 != None: 
       rx_depthN = np.append(0, rx_depthN)
    bg_depth = np.append(0, bg_depth)
    if N7 != None: 
       bg_depthN = np.append(0, bg_depthN)
    dc_depth = np.append(0, dc_depth)
    if N7 != None: 
       dc_depthN = np.append(0, dc_depthN)
    rx_simple_depth = np.append(0, rx_simple_depth)
    if N7 != None: 
       rx_simple_depthN = np.append(0, rx_simple_depthN)
    bg_simple_depth = np.append(0, bg_simple_depth)
    if N7 != None: 
       bg_simple_depthN = np.append(0, bg_simple_depthN)
    dc_simple_depth = np.append(0, dc_simple_depth)
    if N7 != None: 
       dc_simple_depthN = np.append(0, dc_simple_depthN)
    rx_rates = np.append(0, rx_rates)
    if N7 != None: 
       rx_ratesN = np.append(0, rx_ratesN)
    bg_rates = np.append(0, bg_rates)
    if N7 != None: 
       bg_ratesN = np.append(0, bg_ratesN)
    dc_rates = np.append(0, dc_rates)
    if N7 != None: 
       dc_ratesN = np.append(0, dc_ratesN)
    rx_err = np.append(0, rx_err)
    if N7 != None: 
       rx_errN = np.append(0, rx_errN)
    bg_err = np.append(0, bg_err)
    if N7 != None: 
       bg_errN = np.append(0, bg_errN)
    dc_err = np.append(0, dc_err)
    if N7 != None: 
       dc_errN = np.append(0, dc_errN)



    orange_thresh = 0
    red_thresh = 0  
    orange_threshN = 0
    red_threshN = 0  
    if reactivity is not None:

        if dms:
            orange_thresh = 0.2
            red_thresh = 0.4
            if N7 != None:
               orange_threshN = 1.6
               red_threshN = 2.3
        else:
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

        if N7 != None:
           gray_valsN = []
           gray_numsN = []
           gray_errsN = []
           black_valsN = []
           black_numsN = []
           black_errsN = []
           orange_valsN = []
           orange_numsN = []
           orange_errsN = []
           red_valsN = []
           red_numsN = []
           red_errsN = []





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


        if N7 != None:
            
           # If empty, trigger QC
           if (None in reactivityN) or (len(reactivityN) <= 2):
              reactivityN = [0]
              qc_pass_n7="empty"
           # Else, change small reactivities to allow for visualization
           else:
               print("reactivityN: ", reactivityN)
               reactivityN[reactivityN < 0] = .2 

           for i in range(len(reactivityN)):
               if isnan(reactivityN[i]):
                   gray_valsN.append(-1)
                   gray_numsN.append(numN[i])
                   gray_errsN.append(0)
               elif reactivityN[i] < orange_threshN:
                   black_valsN.append(reactivityN[i])
                   black_numsN.append(numN[i])
                   black_errsN.append(stderrN[i])
               elif reactivityN[i] < red_threshN:
                   orange_valsN.append(reactivityN[i])
                   orange_numsN.append(numN[i])
                   orange_errsN.append(stderrN[i])
               else:
                   red_valsN.append(reactivityN[i])
                   red_numsN.append(numN[i])
                   red_errsN.append(stderrN[i])


    
    if dms:
        yMin, ymax = (-0.1, 1.5)
        if N7 != None: 
           yMinN, ymaxN = (-.22, 3.3)
    else:
        yMin, ymax = (-0.5, 4)
    
    left_inches = 0.9
    right_inches = 0.4
    sp_width = len(num)*0.032
    fig_width = max(7,sp_width+left_inches+right_inches)
    fig = plt.figure(figsize=(fig_width,8))

    left_percent = left_inches/fig_width
    right_percent = 1-right_inches/fig_width
    #3 rows, 1 column, which figure it is
    if N7 != None:
       ax1 = plt.subplot(511)
       ax2 = plt.subplot(515)
       ax3 = plt.subplot(512)
       ax4 = plt.subplot(513)
       ax5 = plt.subplot(514)
       plt.subplots_adjust(hspace=0.7, left=left_percent,right=right_percent,top=0.94)
    else:
       ax1 = plt.subplot(311)
       ax2 = plt.subplot(313)
       ax3 = plt.subplot(312)    
       plt.subplots_adjust(hspace=0.5, left=left_percent,right=right_percent,top=0.94)

    near_black = (0,0,1/255.0)


    if reactivity is not  None:
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
    
    mute = False
    if message != None:
       warnings = message.split(",")
       if  "low N7" in warnings:
         mute = True

    if N7 !=  None:
       if reactivity is not None and not mute:
           if len(gray_numsN)>0:
               ax4.bar(gray_numsN,gray_valsN,
                        align="center",
                        width=1.05, color="0.80", edgecolor="0.80",linewidth=0.0)
           if len(black_numsN)>0:
               ax4.bar(black_numsN,black_valsN,
                        align="center",
                        width=1.05, color="black", edgecolor="black",linewidth=0.0,
                        ecolor=near_black,capsize=1)
           if len(orange_numsN)>0:
               ax4.bar(orange_numsN,orange_valsN,
                        align="center",
                        width=1.05, color="pink",edgecolor="pink",linewidth=0.0,
                        ecolor=near_black,capsize=1)
           if len(red_numsN)>0:
               ax4.bar(red_numsN,red_valsN,
                        align="center",
                        width=1.05,color="purple",edgecolor="purple",linewidth=0.0,
                        ecolor=near_black,capsize=1)

    if N7 != None:
         reactivity_graph(ax1, name, yMin, ymax, num, no_mapped, reactivity, orange_thresh, red_thresh, seq, dms = dms, qc_pass = qc_pass, fileout=fileout, N7_reformat = N7_reformat)
         reactivity_graph(ax4, name, yMinN, ymaxN, numN, no_mapped, reactivityN, orange_threshN, red_threshN, seqN, dms = dms, N7 = True, qc_pass_n7=qc_pass_n7, fileout=fileout, message = message, N7_reformat = N7_reformat)
    else: 
         reactivity_graph(ax1, name, yMin, ymax, num, no_mapped, reactivity, orange_thresh, red_thresh, seq, dms = dms, qc_pass = qc_pass, fileout=fileout)

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
    leg = ax2.legend(handles, legend_labels, loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.60)
    xmin, xmax, ymin, ymax = ax2.axis()
    ax2.set_ylim(0, ymax)
    list_ticks = [int(y) for y in ax2.get_yticks()]

    yticks =  np.linspace(0, list_ticks[-1], num = 6) 
    yticks = [int(ytick) for ytick in yticks]

    ax2.set_yticks(yticks) 

    formatted_ticks = []
    for val in yticks:
        formatted_ticks.append(metric_abbreviate(val))
    ax2.set_yticklabels(formatted_ticks)
    if N7 != None:
      ax2.tick_params(axis='y',which='major', labelsize = 9)

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

    if N7:
       if N7_reformat:
         ax2ylabel = ax2.set_ylabel("Read depth", horizontalalignment="left", fontsize=10)
         x, y = ax2ylabel.get_position()
         ax2ylabel.set_position((x, 0))
       else:
         ax2ylabel = ax2.set_ylabel("Read depth", horizontalalignment="left", fontsize=14)
         x, y = ax2ylabel.get_position()
         ax2ylabel.set_position((x, 0))
    else:
      ax2ylabel = ax2.set_ylabel("Read depth", horizontalalignment="left", fontsize=14)
      x, y = ax2ylabel.get_position()
      ax2ylabel.set_position((x, 0))


    ax2.set_xlabel("Nucleotide", fontsize = 12)
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

    if N7 != None:
       good_indicesN = []
       for i in range(len(bg_ratesN)):
           if bg_ratesN[i]<=0.05 or isnan(bg_ratesN[i]):
               good_indicesN.append(i)
       temp_ratesN = [rx_ratesN[i] for i in good_indicesN]
       near_top_rateN = percentile(temp_ratesN,98.0)
       maxesN = [0.32,0.16,0.08,0.04,0.02,0.01]
       ymaxN = maxesN[0]
       for i in range(len(maxesN)):
           if near_top_rateN<maxesN[i]:
               ymaxN = maxesN[i]




    # TODO: put this all into one conditional
    rx_upper = rx_rates + rx_err
    if N7 != None:
       rx_upperN = rx_ratesN + rx_errN
    rx_lower = rx_rates - rx_err
    if N7 != None:
       rx_lowerN = rx_ratesN - rx_errN
    bg_upper = bg_rates + bg_err
    if N7 != None:
       bg_upperN = bg_ratesN + bg_errN
    bg_lower = bg_rates - bg_err
    if N7 != None:
       bg_lowerN = bg_ratesN - bg_errN
    dc_upper = dc_rates + dc_err
    if N7 != None:
       dc_upperN = dc_ratesN + dc_errN
    dc_lower = dc_rates - dc_err
    if N7 != None:
       dc_lowerN = dc_ratesN - dc_errN



    if N7 != None:
       mut_graph(ax3, dms, num, rx_rates, bg_rates, dc_rates, rx_lower, rx_upper, bg_lower, bg_upper, dc_lower, dc_upper, legend_labels, ymax, N7_reformat=N7_reformat)
       mut_graph(ax5, dms, numN, rx_ratesN, bg_ratesN, dc_ratesN, rx_lowerN, rx_upperN, bg_lowerN, bg_upperN, dc_lowerN, dc_upperN, legend_labels, ymaxN, N7 = True, N7_reformat = N7_reformat, message=message)
    else:
       mut_graph(ax3, dms, num, rx_rates, bg_rates, dc_rates, rx_lower, rx_upper, bg_lower, bg_upper, dc_lower, dc_upper, legend_labels, ymax)

    fig.align_ylabels()
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
    ax.text(x1, y, name, transform=ax.transAxes, fontsize=10, horizontalalignment="left", clip_on=False)
    ax.text(x2, y, txt, transform=ax.transAxes, fontsize=10, horizontalalignment="right", clip_on=False)


def write_histograms(shape, stderr,
                     min_depth, max_bg,
                     seq,
                     rx_rate, bg_rate, dc_rate,
                     rx_depth, bg_depth, dc_depth,
                     fileout, name, qc_pass):


    # limit rate histograms to uppercase sequence
    ui = np.array([False]*len(seq))
    for i, c in enumerate(seq):
        if not c.islower():
            ui[i] = True
    if rx_rate is not None:
        rx_rate = np.array(rx_rate)[ui]
        rx_depth = np.array(rx_depth)[ui]
    if bg_rate is not None:
        bg_rate = np.array(bg_rate)[ui]
        bg_depth = np.array(bg_depth)[ui]
    if dc_rate is not None:
        dc_rate = np.array(dc_rate)[ui]
        dc_depth = np.array(dc_depth)[ui]

    # limit rate histograms to positions included in final reactivity profile
    if shape is not None: # If shape is None this functionality will break
        fin_mask = np.isfinite(shape)
        if sum(ui) > 0:
            fin_mask = np.isfinite(shape)[ui]
        if rx_rate is not None:
            rx_rate = rx_rate[fin_mask]
            rx_depth = rx_depth[fin_mask]
        if bg_rate is not None:
            bg_rate = bg_rate[fin_mask]
            bg_depth = bg_depth[fin_mask]
        if dc_rate is not None:
            dc_rate = dc_rate[fin_mask]
            dc_depth = dc_depth[fin_mask]

    tlabel_size = 10

    legend_labels = []
    if rx_rate is not None:
        legend_labels.append("Modified")
    if bg_rate is not None:
        legend_labels.append("Untreated")
    if dc_rate is not None:
        legend_labels.append("Denatured")

    fig, axes = plt.subplots(nrows=2, ncols=3)
    ax1, ax2, ax3 = axes[0,:]
    ax4 = axes[1,2]
    fig.set_size_inches(10,6)
    fig.subplots_adjust(bottom=0.1, top=0.85, wspace=0.5, hspace=0.75, left=0.08, right=0.95)

    for ax in axes[1,0:2]:
        ax.set_visible(False)

    if len(name)<20:
        title = plt.suptitle(name,fontsize=18,horizontalalignment="left", x=0.02)
    if not qc_pass:
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

        if return_flag:
            plt.savefig(fileout)
            return
 
 
    num_bins = 30
    # use fewer bins for short RNAs
    if len(rx_depth)<500:
        num_bins = 10
    int_bins = range(0,num_bins+1,1)
    

    ptiles = []
    if rx_rate is not None and len(rx_rate)>10:
        ptiles.append(percentile(rx_rate, 90.0))
    if bg_rate is not None and len(bg_rate)>10:
        ptiles.append(percentile(bg_rate, 90.0))
    if dc_rate is not None and len(dc_rate)>10:
        ptiles.append(percentile(dc_rate, 90.0))
    if len(ptiles) == 0:
        max90percentile = 0
    else:
        max90percentile = max(ptiles)

    # adjust rate axis for higher mutation rates (e.g. from DMS or other highly mutagenic reagents)
    bin_max = 0.03
    if max90percentile < 0.01:
        bin_max = 0.008
    bin_width = bin_max/(num_bins)
    float_bins = [b*bin_width for b in int_bins]

    ax1.set_ylabel("Normalized\nnucleotide count",fontsize=13)
    ax1.set_xlabel("Mutation rate (%)", fontsize=13)
    ax1.set_title('Mutation rates', x=0.5,y=1.08)

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
        draw_percentile(ax1, rx_rate, 95.0, "black", x1=-0.2, x2=1, y=-1, percentage=True, name="Modified sample:\n 95th percentile rate:")
        draw_percentile(ax1, rx_rate, 50.0, "black", x1=-0.2, x2=1, y=-1.1, percentage=True, name=" Median rate:")

    if "Untreated" in legend_labels and len(bg_rate)>10:
        draw_percentile(ax1, bg_rate, 95.0, "black", x1=-0.2, x2=1, y=-1.4, percentage=True, name="Untreated sample:\n 95th percentile rate:")
        draw_percentile(ax1, bg_rate, 50.0, "black", x1=-0.2, x2=1, y=-1.5, percentage=True, name=" Median rate:")

    if "Denatured" in legend_labels and len(dc_rate)>10:
        draw_percentile(ax1, dc_rate, 95.0, "black", x1=-0.2, x2=1, y=-1.8, percentage=True, name="Denatured sample:\n 95th percentile rate:")
        draw_percentile(ax1, dc_rate, 50.0, "black", x1=-0.2, x2=1, y=-1.9, percentage=True, name=" Median rate:")

    # Plot ln(RXrate/BGrate) histogram in lower right panel
    if "Modified" in legend_labels and "Untreated" in legend_labels:
        vals = np.log(rx_rate/bg_rate)
        vals = vals[np.isfinite(vals)]
        # typical experiments range between -1 and 7
        pdf, bins, patches = ax4.hist(vals, bins=num_bins, range=(-1,7), histtype='step', lw=1.5,
                                      color='black')
        max_y = max([x[1] for x in patches[0].get_xy()])
        # rescale stepplot to max out at 1
        patches[0].set_xy([[c[0],c[1]/max_y] for c in patches[0].get_xy()])
        ax4.set_xlabel("ln(modified rate / untreated rate)", fontsize=13)
        ax4.set_ylabel("Normalized\nnucleotide count", fontsize=13)
        ax4.set_title("Raw background-corrected\nrate distribution", x=0.5, y=1.04)
        ax4.spines["right"].set_visible(False)
        ax4.spines["top"].set_visible(False)
        ax4.tick_params(axis='both', labelsize=tlabel_size)
        ax4.get_xaxis().tick_bottom()
        ax4.get_yaxis().tick_left()
        ax4.set_ylim([0, 1])
        ax4.set_xticks(list(range(-1,8)))


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

    ymax = 1
    ax3.set_ylim([0,ymax])

    ax3.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax3.get_yaxis().tick_left()

    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    for line in ax3.get_yticklines() + ax3.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(1)

    plt.savefig(fileout)

 
def write_histograms_dms(dms, stderr,
                         min_depth, max_bg,
                         seq,
                         rx_rate, bg_rate, dc_rate,
                         rx_depth, bg_depth, dc_depth,
                         fileout, name, qc_pass, N7file=None, dms_N7 = None,
                         rx_rate_N7=None, bg_rate_N7=None, dc_rate_N7=None,
                         rx_depth_N7=None, bg_depth_N7=None, dc_depth_N7=None):
    
    
    seq = np.array(list(seq))

    tlabel_size = 10

    legend_labels = []
    if rx_rate is not None:
        legend_labels.append("Modified")
    if bg_rate is not None:
        legend_labels.append("Untreated")
    if dc_rate is not None:
        legend_labels.append("Denatured")
    
    if N7file:
        if rx_rate is not None and bg_rate is not None:
            mkdiffplot=True
            fig, axes = plt.subplots(nrows=5, ncols=4)
            fig.set_size_inches(10,12.5)
        else:
            mkdiffplot=False
            fig, axes = plt.subplots(nrows=5, ncols=3)
            fig.set_size_inches(7.5,12.5)

    else:
        if rx_rate is not None and bg_rate is not None:
            mkdiffplot=True
            fig, axes = plt.subplots(nrows=4, ncols=4)
            fig.set_size_inches(10,10)
        else:
            mkdiffplot=False
            fig, axes = plt.subplots(nrows=4, ncols=3)
            fig.set_size_inches(7.5,10)


    fig.subplots_adjust(bottom=0.1, top=0.85, wspace=0.25, hspace=0.45, left=0.08, right=0.95)

    if len(name)<20:
        title = plt.suptitle(name,fontsize=18,horizontalalignment="left", x=0.02)
    if not qc_pass:
        if len(name)<20:
            offset = len(name)*12
        else:
            offset = 0
        msg = "Note: possible data quality issue - see log file"
        return_flag = False        
        if bg_rate is None and rx_rate is None and dc_rate is None:
            return_flag = True
            msg = "ERROR: no reads mapped to this RNA"
        txt = plt.text(0.1, 0.95,
                 msg,
                 ha='left', va='top',
                 fontsize=16, color='red',
                 transform=fig.transFigure)

        if return_flag:
            plt.savefig(fileout)
            return


    def _plothist_(ax, mask, rx, bg, dc, nbins, GN7 = False):

        minval, maxval = 1, 0
        
        for data in (rx, bg, dc):
            
            if data is None:
                continue

            m = min(data[mask])
            if m < minval:
                minval=m
            m = 1.2*np.percentile(data[mask], 95)
            if m>maxval:
                maxval=m
            
        bins = np.linspace(minval,maxval,nbins)
        inc = 0



        if GN7:
            if rx is not None:
                ax.hist(rx[mask], bins=bins, histtype='step',lw=1.5,color="purple", zorder=4)
                if max(rx[mask]) > maxval:
                    plt.text(1.0, 0.3, 'max={:.3f} not shown'.format(max(rx[mask])), transform=ax.transAxes, fontsize=6, color="purple", ha='right', style='italic')
                    inc+=0.1
            if bg is not None and sum(bg>0) > 0:
                ax.hist(bg[mask], bins=bins, histtype='step',lw=1.5,color="pink", zorder=3)
                if max(bg[mask]) > maxval:
                    plt.text(1.0, 0.3-inc, 'max={:.3f} not shown'.format(max(bg[mask])), transform=ax.transAxes, fontsize=6, color="pink", ha='right', style='italic')
                    inc+=0.1
            if dc is not None and sum(dc>0) > 0:
                ax.hist(dc[mask], bins=bins, histtype='step',lw=1.5,color=dc_color, zorder=2)
                if max(dc[mask]) > maxval:
                    plt.text(1.0, 0.3-inc, 'max={:.3f} not shown'.format(max(dc[mask])), transform=ax.transAxes, fontsize=6, color=dc_color, ha='right', style='italic')

        else:
            if rx is not None:
                ax.hist(rx[mask], bins=bins, histtype='step',lw=1.5,color=rx_color, zorder=4)
                if max(rx[mask]) > maxval:
                    plt.text(1.0, 0.3, 'max={:.3f} not shown'.format(max(rx[mask])), transform=ax.transAxes, fontsize=6, color=rx_color, ha='right', style='italic')
                    inc+=0.1
            if bg is not None and sum(bg>0) > 0:
                ax.hist(bg[mask], bins=bins, histtype='step',lw=1.5,color=bg_color, zorder=3)
                if max(bg[mask]) > maxval:
                    plt.text(1.0, 0.3-inc, 'max={:.3f} not shown'.format(max(bg[mask])), transform=ax.transAxes, fontsize=6, color=bg_color, ha='right', style='italic')
                    inc+=0.1
            if dc is not None and sum(dc>0) > 0:
                ax.hist(dc[mask], bins=bins, histtype='step',lw=1.5,color=dc_color, zorder=2)
                if max(dc[mask]) > maxval:
                    plt.text(1.0, 0.3-inc, 'max={:.3f} not shown'.format(max(dc[mask])), transform=ax.transAxes, fontsize=6, color=dc_color, ha='right', style='italic')


    if N7file:
        NT_iterable = ['A','C','U','G','GN7']
    else:
        NT_iterable = ['A','C','U','G']

    for ntindex, nt in enumerate(NT_iterable):


        # In the future this logic should be incorporated in a more integrated way
        # as this is a bit hacky
        if nt == "GN7":
            rx_rate = rx_rate_N7
            bg_rate = bg_rate_N7
            dc_rate = dc_rate_N7
            dms = dms_N7

        if nt == "GN7":
            mask = (seq == 'G')
        else: 
            mask = (seq == nt)

        for rate in (rx_rate, bg_rate, dc_rate):
            if rate is not None:
                mask = mask & np.isfinite(rate)
        
        if mkdiffplot:
            ax1, ax2, ax3, ax4 = axes[ntindex,:]
        else:
            ax1, ax3, ax4 = axes[ntindex,:]


        plt.text(-0.3,1.0, nt, transform=ax1.transAxes, fontsize=20, weight='bold')

        ax1.get_xaxis().tick_bottom() 
        ax1.get_yaxis().tick_left()

        ax1.spines["right"].set_visible(False)
        ax1.spines["top"].set_visible(False)
        for line in ax1.get_yticklines() + ax1.get_xticklines():
            line.set_markersize(7)
            line.set_markeredgewidth(1)

        ax1.tick_params(axis='both', labelsize=tlabel_size)

        if sum(mask) == 0: # No NT rates to interact with / plot.
                           # Skip to avoid error
            continue

        if nt == "GN7":
            _plothist_(ax1, mask, rx_rate, bg_rate, dc_rate, 20, True)
        else:
            _plothist_(ax1, mask, rx_rate, bg_rate, dc_rate, 20)

        ax1.set_xticklabels(['{:.3g}'.format(x) for x in ax1.get_xticks()], rotation = 90)
        
        if nt=='A' or nt=='GN7':
            leg = ax1.legend(legend_labels, framealpha=0.75, fontsize=8, markerscale=0.5)
            for l in leg.get_lines():
                l.set_linewidth(1.5)
                l.set_alpha(1.0)

        # plot the reactivity difference
        if mkdiffplot:
            diff = rx_rate-bg_rate
            p95 = 1.2*np.percentile(diff[mask], 95)
            if p95 < 0:
                # In cases where all rates are negative, ie high background, low number NTs, or both
                # p95 will be transformed to be 20% smaller than the 95th percentile.
                # So instead we need to make it 20% larger. IE move it 20% in the positive direction
                # instead of negative to avoid error due to the range having a max smaller than the min.
                p95 = .8*np.percentile(diff[mask], 95)
                

            ax2.hist(diff[mask], bins=20, range=(min(diff[mask]),p95), histtype='step', lw=1.5, color='black')
            ax2.spines["right"].set_visible(False)
            ax2.spines["top"].set_visible(False)
            ax2.tick_params(axis='both', labelsize=tlabel_size)
            ax2.get_xaxis().tick_bottom()
            ax2.get_yaxis().tick_left()

            ax2.set_xticklabels(['{:.3g}'.format(x) for x in ax2.get_xticks()], rotation = 90)
            
            if max(diff[mask]) > p95:
                plt.text(1.0, 0.3, 'max={:.3f} not shown'.format(max(diff[mask])), transform=ax2.transAxes, fontsize=6, color='k', ha='right', style='italic')


        if nt == "GN7":
            xticks = np.linspace(0.0,3.3,11)
        else:
            xticks = np.linspace(-0.4,1.6,11)
        xlabels = ["{:.0f}".format(val) if val==int(val) else "{:.1f}".format(val) for val in xticks]
        
        ax3.tick_params(axis='both', labelsize=tlabel_size)
        ax3.set_xticks(xticks)
        ax3.set_xticklabels(xlabels, rotation=90)

        if dms is not None:
            finmask = mask & np.isfinite(dms)
            if sum(finmask) > 0:

                if nt == "GN7":
                    
                    
                    d = ax3.hist(dms[finmask], bins=20, range=(0, 3.3), histtype='step', lw=1.5, color="black", zorder=4)

                    #proceeds_purine_mask = [nt for indx, nt_elem in enumerate(seq)]
                    # If num purine or pyrimidine too low NFAC is pooled
                    g_pur = ( [seq[nt_index + 1] in ["A", "G"] for nt_index, nt_value in enumerate(seq[:-1])] + [False] ) & finmask
                    g_pyr = ( [seq[nt_index + 1] in ["C", "T", "U"] for nt_index, nt_value in enumerate(seq[:-1])] + [False] ) & finmask

                    exp2_dms = np.exp2(dms) # Reverse log2 transformation so norm fac can be calculated
                    if sum(g_pur) < 15 or sum(g_pyr) < 15:
                        idx = np.argmax(exp2_dms[finmask])
                        if mkdiffplot:
                            nfac = diff[finmask][idx]*exp2_dms[finmask][idx]
                            ax2.axvline(x=nfac, ls='--', lw=1.5, color='0.5')
                            plt.text(1.0,1.0, 'NormFac={:.3g}'.format(nfac), transform=ax2.transAxes, fontsize=8, color='0.5', ha='right')
                        else:
                            nfac = rx_rate[finmask][idx]*exp2_dms[finmask][idx]
                            ax1.axvline(x=nfac, ls='--', lw=1.5, color='0.5')
                            plt.text(1.0,1.0, 'NormFac={:.3g}'.format(nfac), transform=ax1.transAxes, fontsize=8, color='0.5', ha='right')
                    # Else, seperate NFAC for G_PUR and G_PYR
                    else:
                        g_pur_idx = np.argmax(exp2_dms[g_pur])
                        g_pyr_idx = np.argmax(exp2_dms[g_pyr])
                        if mkdiffplot:
                            g_pur_nfac = exp2_dms[g_pur][g_pur_idx] * diff[g_pur][g_pur_idx]
                            g_pyr_nfac = exp2_dms[g_pyr][g_pyr_idx] * diff[g_pyr][g_pyr_idx]
                            ax2.axvline(x = g_pur_nfac, ls='--', lw=1.5, color='0.5')
                            ax2.axvline(x = g_pyr_nfac, ls='--', lw=1.5, color='0.5')
                            plt.text(1.0,1.0, 'GPUR={:.3g}, GPYR={:.3g}'.format(g_pur_nfac, g_pyr_nfac), transform=ax2.transAxes, fontsize=8, color='0.5', ha='right')

                        else:
                            g_pur_nfac = exp2_dms[g_pur][g_pur_idx] * rx_rate[g_pur][g_pur_idx]
                            g_pyr_nfac = exp2_dms[g_pyr][g_pyr_idx] * rx_rate[g_pyr][g_pyr_idx]
                            ax1.axvline(x = g_pur_nfac, ls='--', lw=1.5, color='0.5')
                            ax1.axvline(x = g_pyr_nfac, ls='--', lw=1.5, color='0.5')
                            plt.text(1.0,1.0, 'GPUR={:.3g}, GPYR={:.3g}'.format(g_pur_nfac, g_pyr_nfac), transform=ax1.transAxes, fontsize=8, color='0.5', ha='right')

                else: 
                    d = ax3.hist(dms[finmask], bins=20, range=(-0.4,1.6), histtype='step', lw=1.5, color="black", zorder=4)
                    #plot normalization factor
                    idx = np.argmax(dms[finmask])
                    if mkdiffplot:
                        nfac = diff[finmask][idx]/dms[finmask][idx]
                        ax2.axvline(x=nfac, ls='--', lw=1.5, color='0.5')
                        plt.text(1.0,1.0, 'NormFac={:.3g}'.format(nfac), transform=ax2.transAxes, fontsize=8, color='0.5', ha='right')
                    else:
                        nfac = rx_rate[finmask][idx]/dms[finmask][idx]
                        ax1.axvline(x=nfac, ls='--', lw=1.5, color='0.5')
                        plt.text(1.0,1.0, 'NormFac={:.3g}'.format(nfac), transform=ax1.transAxes, fontsize=8, color='0.5', ha='right')


        ax3.get_xaxis().tick_bottom()   # remove unneeded ticks
        ax3.get_yaxis().tick_left()

        ax3.spines["right"].set_visible(False)
        ax3.spines["top"].set_visible(False)
        for line in ax3.get_yticklines() + ax3.get_xticklines():
            line.set_markersize(7)
            line.set_markeredgewidth(1)


        if nt == "GN7":
            _plothist_(ax4, mask, rx_depth, bg_depth, dc_depth, 20, True)
        else:
            _plothist_(ax4, mask, rx_depth, bg_depth, dc_depth, 20)

        ax4.get_xaxis().tick_bottom()   # remove unneeded ticks
        ax4.get_yaxis().tick_left()
        ax4.tick_params(axis='both', labelsize=tlabel_size)

        xticks = [int(x) for x in ax4.get_xticks()]
        formatted_ticks = []
        for val in xticks:
            formatted_ticks.append(metric_abbreviate(val))
        ax4.set_xticklabels(formatted_ticks, rotation=90)

        ax4.spines["right"].set_visible(False)
        ax4.spines["top"].set_visible(False)

        for line in ax4.get_yticklines() + ax4.get_xticklines():
            line.set_markersize(7)
            line.set_markeredgewidth(1)

        
        
        ax1.set_ylabel("Nucleotide count",fontsize=11)
        
        if ntindex==0:
            ax1.set_title('Mutation rates', x=0.5,y=1.08)
            ax3.set_title("Normalized\nreactivity distribution", x=0.5,y=1.08)
            ax4.set_title("Read depths", x=0.5,y=1.08)
            if mkdiffplot:
                ax2.set_title("Background corrected\nrate distribution", x=0.5, y=1.08)

        if (ntindex==3 and not N7file) or (ntindex==4 and N7file) :
            ax1.set_xlabel("Mutation rate", fontsize=11)
            ax3.set_xlabel("Normalized reactivity", fontsize=11)
            ax4.set_xlabel("Effective read depth", fontsize=11)
            if mkdiffplot:
                ax2.set_xlabel("Modified - Untreated Rate", fontsize=11)


    plt.savefig(fileout)



def load_tab(filename):
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

    for i in range(length):
        line = f.readline()
        s = line.strip().split('\t')
        for j in range(len(headers)):
            try:
                data[headers[j]][i] = eval(remapped_types[j])(s[j])
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

    h = "Tab-delimited file summarizing N7 reactivity"
    h += " profile, mutation rates, and sequencing depth. "
    parser.add_argument("--N7file", help=h, required=False, type=str)

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

    h = "Amplicon primer pair sequences and locations (to exclude from mutation rate histogram plots)"
    parser.add_argument("--primers", help=h, type=str)
    
    h = "Reactivities are DMS rather than default SHAPE"
    parser.add_argument("--dms", action="store_true", default=False, help=h)
    
    p = parser.parse_args(sys.argv[1:])

    primers = []
    if p.primers is not None and p.primers != "":
        primers = load_primers(p.primers)

    d = load_tab(p.infile)
    if(p.N7file != None):
       dn = load_tab(p.N7file)

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
        if(p.N7file != None):
           profileN = dn["Norm_profile"]
           stderrN = dn["Norm_stderr"] 
    else:
        profile = None
        stderr = None
        if(p.N7file != None):
           profileN = None
           stderrN = None 

    if p.title is not None:
        title = p.title
        if(p.N7file != None):
           titleN = p.title
        
    else:
        title = ""
        if(p.N7file != None):
           titleN = ""

    # mask out amplicon primer pair site sequences before running QC checks
    masked_sequence = d["Sequence"].copy()
    if(p.N7file != None):
       masked_sequenceN = dn["Sequence"].copy()

    masked = np.zeros(masked_sequence.size, dtype=bool)
    if(p.N7file != None):
        maskedN = np.zeros(masked_sequenceN.size, dtype=bool)

    for pair in primers:
        for primer in [pair.fw, pair.rv]:
            for i in range(primer.left-1, primer.right):
                try:
                    masked[i] = True
                    if(p.N7file != None):
                       masked[i] = True
                except IndexError:
                    pass
    for i in range(len(masked)):
        if masked[i]:
            masked_sequence[i] = masked_sequence[i].lower()
            if(p.N7file != None):
               masked_sequenceN[i] = masked_sequenceN[i].lower()
    
    
    if p.dms:
        if "Norm_profile" in d:
            norm_p = d["Norm_profile"]
        else:
            norm_p = None

        qc_pass = qc_stats_dms(masked_sequence,
                               d["Modified_rate"],
                               d["Untreated_rate"],
                               d["Denatured_rate"],
                               d["Modified_effective_depth"],
                               d["Untreated_effective_depth"],
                               d["Denatured_effective_depth"],
                               norm_p,
                               p.mindepth, p.maxbg,
                               p.min_depth_pass,
                               p.max_high_bg,
                               p.min_positive,
                               p.min_high_mut)
    #TODO: move N7 QC here.

    else:
        qc_pass = qc_stats(masked_sequence,
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
         
        if(p.N7file == None):
            render_profiles(
               d["Nucleotide"], d["Sequence"], profile, stderr,
               d["Modified_rate"], d["Untreated_rate"], d["Denatured_rate"],
               d["Modified_effective_depth"], d["Untreated_effective_depth"], d["Denatured_effective_depth"],
               d["Modified_read_depth"], d["Untreated_read_depth"], d["Denatured_read_depth"],
               p.plot, title, qc_pass, dms=p.dms)

        else:
            #Logic to scrape / N7 quality control message produced in norm profile
            n7_split = p.N7file.split("/")
            prefix = "/".join(n7_split[0:-1])
            suffix = n7_split[-1]
            n7_message_file = prefix + "/." + suffix + "_n7message.txt"

            try:
               n7_inp_file = open(n7_message_file, "r")
               lines = [line for line in n7_inp_file]
               n7_inp_file.close()
               message = lines[0]

            except:
               message = None


            render_profiles(
               d["Nucleotide"], d["Sequence"], profile, stderr,
               d["Modified_rate"], d["Untreated_rate"], d["Denatured_rate"],
               d["Modified_effective_depth"], d["Untreated_effective_depth"], d["Denatured_effective_depth"],
               d["Modified_read_depth"], d["Untreated_read_depth"], d["Denatured_read_depth"],
               p.plot, title, qc_pass, 
               dn["Nucleotide"], dn["Sequence"], profileN, stderrN,
               dn["Modified_rate"], dn["Untreated_rate"], dn["Denatured_rate"],
               dn["Modified_effective_depth"], dn["Untreated_effective_depth"], dn["Denatured_effective_depth"],
               dn["Modified_read_depth"], dn["Untreated_read_depth"], dn["Denatured_read_depth"],
               dms=p.dms, N7=p.N7file, message=message 
               )



    if p.hist is not None:
        if p.N7file:
            write_histograms_dms(
                profile, stderr,
                p.mindepth, p.maxbg,
                masked_sequence,
                d["Modified_rate"], d["Untreated_rate"], d["Denatured_rate"],
                d["Modified_effective_depth"], d["Untreated_effective_depth"], d["Denatured_effective_depth"],
                p.hist, title, qc_pass, p.N7file,
                profileN, dn["Modified_rate"], dn["Untreated_rate"], dn["Denatured_rate"],
                dn["Modified_effective_depth"], dn["Untreated_effective_depth"], dn["Denatured_effective_depth"])
        elif p.dms:
            write_histograms_dms(
                profile, stderr,
                p.mindepth, p.maxbg,
                masked_sequence,
                d["Modified_rate"], d["Untreated_rate"], d["Denatured_rate"],
                d["Modified_effective_depth"], d["Untreated_effective_depth"], d["Denatured_effective_depth"],
                p.hist, title, qc_pass)
        else:
            write_histograms(
                profile, stderr,
                p.mindepth, p.maxbg,
                masked_sequence,
                d["Modified_rate"], d["Untreated_rate"], d["Denatured_rate"],
                d["Modified_effective_depth"], d["Untreated_effective_depth"], d["Denatured_effective_depth"],
                p.hist, title, qc_pass)
