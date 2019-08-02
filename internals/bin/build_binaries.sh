#!/usr/bin/env bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"
cd ${BASE_DIR}

echo "Building ShapeMapper executables..."
#rm -rf build
mkdir -p build
cd build
cmake ..
make

source ${BASE_DIR}/internals/paths/bin_paths.sh

if [ -z $(which shapemapper_read_trimmer) ] || \
   [ -z $(which shapemapper_mutation_parser) ] || \
   [ -z $(which shapemapper_mutation_counter) ]; then
    msg="Error building ShapeMapper executables."
    echo "$msg"
    exit 1
else
    msg="ShapeMapper executables successfully built."
    echo "$msg"
fi

