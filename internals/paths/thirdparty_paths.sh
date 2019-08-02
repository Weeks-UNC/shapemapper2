#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"
CONDA_PATH=${BASE_DIR}/internals/thirdparty/miniconda


# python, bbmerge, graphviz (dot), STAR
export PATH="${CONDA_PATH}/bin:${PATH}"

# perl
export PERL5LIB="${CONDA_PATH}/lib/perl5/5.22.0"
export PERL5LIB="${CONDA_PATH}/lib/perl5/5.22.0/x86_64-linux-thread-multi:${PERL5LIB}"
export PATH="${CONDA_PATH}:${PATH}"

# bowtie2
export PATH="${BASE_DIR}/internals/thirdparty/bowtie2:${PATH}"

# pipeviewer
export PATH="${BASE_DIR}/internals/thirdparty/pipeviewer/bin:${PATH}"

# ghostscript
export PATH="${BASE_DIR}/internals/thirdparty/ghostscript:${PATH}"
