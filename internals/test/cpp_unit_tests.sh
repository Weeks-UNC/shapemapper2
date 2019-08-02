#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

# Run C++ unit tests (Google test framework)

#set -e # exit on first error (if any)
CURRENT_DIR="$(pwd)"

# Find the parent folder of this script,
# resolving (possibly nested) symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$BASE_DIR/$SOURCE"
done
BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && cd ../.. && pwd )"


source ${BASE_DIR}/internals/paths/bin_paths.sh

# disable core dumps
ulimit -c 0

err="ERROR: error running c++ unit test executable.\nFAILURE"

test_read_trimmer "${BASE_DIR}"
if [[ $? != 0 ]]; then
    echo -e "${err}"
    exit $?
fi

test_mutation_parser "${BASE_DIR}"
if [[ $? != 0 ]]; then
    echo -e "${err}"
    exit $?
fi

test_mutation_counter "${BASE_DIR}"
if [[ $? != 0 ]]; then
    echo -e "${err}"
    exit $?
fi

test_histogram
if [[ $? != 0 ]]; then
    echo -e "${err}"
    exit $?
fi
