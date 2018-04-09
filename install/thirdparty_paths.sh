#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
CONDA_PATH=${BASE_DIR}/thirdparty/miniconda


# python, bbmerge, graphviz (dot), STAR
export PATH="${CONDA_PATH}/bin:${PATH}"

# bowtie2
export PATH="${CONDA_PATH}:${PATH}"

# pipeviewer
export PATH="${BASE_DIR}/thirdparty/pipeviewer/bin:${PATH}"
