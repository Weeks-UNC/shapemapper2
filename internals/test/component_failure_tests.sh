#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

# Tests for pipeline error detection at specific modules

#set -e # exit on first error (if any)

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

#cd ${BASE_DIR}/example_data
cd ${BASE_DIR}/internals/test/data


echo '[==========] Fail testing pipeline: 1-sample, correct-seq, parallel'
python3 ${BASE_DIR}/internals/test/fail_tester.py \
--target TPP.fa \
--out fail_tester_results \
--temp fail_tester_temp \
--log fail_tester_results/log.txt \
--name fail_tester \
--overwrite \
--min-depth 10 \
--modified --folder TPPplus \
--correct-seq --folder TPPminus

rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi


echo -e "\n\n"
echo '[==========] Fail testing pipeline: 1-sample, serial'
python3 ${BASE_DIR}/internals/test/fail_tester.py \
--target TPP.fa \
--out fail_tester_results \
--temp fail_tester_temp \
--log fail_tester_results/log.txt \
--name fail_tester_serial \
--overwrite \
--min-depth 10 \
--modified --folder TPPplus \
--serial

rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

echo -e "\n\n"
echo '[==========] Fail testing pipeline: 1-sample, correct-seq, parallel, STAR aligner'
python3 ${BASE_DIR}/internals/test/fail_tester.py \
--target TPP.fa \
--out fail_tester_results \
--temp fail_tester_temp \
--log fail_tester_results/log.txt \
--name fail_tester_star \
--star-aligner \
--overwrite \
--min-depth 10 \
--modified --folder TPPplus \
--correct-seq --folder TPPminus

rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

rm -rf core.* # remove any core dumps
rm -rf fail_tester_results
rm -rf fail_tester_temp
