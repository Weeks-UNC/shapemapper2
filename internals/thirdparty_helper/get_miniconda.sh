#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

# Check for local miniconda, download if necessary

set -e # exit on first error (if any)
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"

cd ${BASE_DIR}/internals

CONDA_PATH=${BASE_DIR}/internals/thirdparty/miniconda
export PATH=${CONDA_PATH}/bin:$PATH
# Check for local miniconda, download and extract if not present
if [ -f "${BASE_DIR}/internals/thirdparty/miniconda/bin/conda" ]; then
    echo "Miniconda already present locally"
else
    # more recent miniconda3 releases have apparently moved to
    # requiring at least GLIBC 2.6
    #p='https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh'
    p='https://repo.continuum.io/miniconda/Miniconda3-4.3.21-Linux-x86_64.sh'
    mkdir -p thirdparty
    cd thirdparty

    echo "Downloading base miniconda environment"
    curl -L "${p}" -o miniconda.sh

    unset PYTHONPATH

    rm -rf "${CONDA_PATH}"
    bash miniconda.sh -b -p "${CONDA_PATH}"
fi

