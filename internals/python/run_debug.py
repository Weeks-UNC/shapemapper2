#!/usr/bin/env python3
"""
Pipeline runner, for easier debugging from PyCharm. Not actually
used in normal ShapeMapper execution.
"""
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #


import os, subprocess
from cli import run

this_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(this_dir, "..")

# load executable paths from shapemapper/install/bin_paths.sh
# - no way to do this directly from a child process, so print and parse
cmd = '. "{base_dir}/install/bin_paths.sh"'.format(base_dir=base_dir)
cmd += ' && echo "${PATH}"'
s = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
s = s.decode("utf-8")
os.environ["PATH"] = s


os.chdir(os.path.join(base_dir,"example_data"))


# FIXME: run with a smaller --max-paired-fragment-length and --render-mutations
#        to generate usable images for docs

args = """shapemapper
--target 16S.fa 23S.fa
--out example_results
--log example_results/log.txt
--name ribosome-for-docs
--render-mutations
--overwrite
--max-paired-fragment-length 200
--modified --folder ribosome_plus
""".split()

args = """shapemapper
--target 16S.fa 23S.fa
--out example_results
--log example_results/log.txt
--name ribosome-for-docs
--render-mutations
--overwrite
--verbose
--max-paired-fragment-length 200
--max-pages 600
--modified --folder ../test/data/ribosome_plus
""".split()



args = """shapemapper
--target ../test/data/TPP.fa
--out example_results
--log example_results/log.txt
--name TPP
--modified --folder TPPplus
--render-flowchart
--amplicon
--overwrite
""".split()


args = """shapemapper
--target ../test/data/TPP_masked_primers.fa
--out example_results
--log example_results/log.txt
--name TPP
--modified --folder TPPplus
--render-flowchart
--amplicon
--star-aligner
--genomeSAindexNbase 5
--rerun-on-star-segfault
--rerun-genomeSAindexNbase 3
--overwrite
""".split()

args = """shapemapper
--target ../test/data/primer_pairs/multiple_pairs/target_regions.fa
--primers ../test/data/primer_pairs/multiple_pairs/primers
--out example_results
--log example_results/log.txt
--name multiple-amplicons
--correct-seq --folder ../test/data/primer_pairs/multiple_pairs/minus
--overwrite
--render-flowchart
""".split()

#--primers ../test/data/primer_pairs/multiple_pairs/primers

# for generating paired, non-overlapping reads
#--min-qual-to-trim 30
#--window-to-trim 1

args = """shapemapper
--target ../test/data/primer_pairs/RMRP/Rmrp.fa
--out example_results
--name RMRP
--overwrite
--min-depth 1000
--max-paired-fragment-length 200
--amplicon
--serial
--verbose
--preserve-order
--render-mutations
--output-parsed
--output-aligned
--max-primer-offset 20
--modified --folder ../test/data/primer_pairs/RMRP/plus
""".split()




#---------------------------------------------------------

args = """shapemapper
--target ../test/data/primer_pairs/multiple_pairs/target_regions.fa
--primers ../test/data/primer_pairs/multiple_pairs/primers
--out /home/wahoosteve/Documents/group_meeting_9-18-18/tiled-amplicons
--temp /home/wahoosteve/Documents/group_meeting_9-18-18/tiled-amplicons
--log /home/wahoosteve/Documents/group_meeting_9-18-18/tiled-amplicons/log.txt
--output-parsed
--render-mutations
--name tiled-amplicons
--modified --folder ../test/data/primer_pairs/multiple_pairs_more_depth/plus
--untreated --folder ../test/data/primer_pairs/multiple_pairs_more_depth/minus
--overwrite
--render-flowchart
""".split()



#-----------------------------------------------------------



# ../test/data/TPP.fa is unmasked (no lowercase)

#--verbose
#--preserve-order
#--render-mutations

#--correct-seq --folder TPPminus
#--modified --folder TPPplus
#--untreated --folder TPPminus
#--denatured --folder TPPdenat
#--ambig

#shapemapper \
#--target TPP.fa \
#--out example_results \
#--log example_results/log.txt \
#--name TPP_three_samples \
#--modified --folder TPPplus #\
##--untreated --folder TPPminus \
##--denatured --folder TPPdenat





