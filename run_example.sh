#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

# Run simple ShapeMapper pipeline on a small subset of example data

set -e # exit on first error (if any)

# Find the parent folder of this script,
# resolving (possibly nested) symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$BASE_DIR/$SOURCE"
done
BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

export PATH=${BASE_DIR}:${PATH}

cd ${BASE_DIR}/example_data

shapemapper \
--name "example-results" \
--target TPP.fa \
--amplicon \
--overwrite \
--min-depth 4000 \
--modified --folder TPPplus \
--untreated --folder TPPminus \
--denatured --folder TPPdenat
