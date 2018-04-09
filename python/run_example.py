"""
Python example runner, for easier debugging from PyCharm. Not actually
used in normal ShapeMapper execution.
"""
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

# FIXME: load paths from shapemapper/install/bin_paths.sh

import os
from cli import run

this_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(this_dir, "..")

os.chdir(os.path.join(base_dir,"example_data"))

#os.environ["PATH"] += os.pathsep + os.path.join(base_dir)

args = """shapemapper
--target 16S.fa 23S.fa
--out example_results
--log example_results/log.txt
--name ribosome
--modified --folder ribosome_plus
""".split()

run(args)

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





