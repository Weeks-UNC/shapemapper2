#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

# Analyze subset of bacterial ribosome dataset and check
# that area under ROC curve metric is above a threshold

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

cd ${BASE_DIR}/internals/test/data

mkdir -p shapemapper_temp


#echo "bowtie2 version:"
#bowtie2 --version
#echo -e '\n'
#echo "${PATH}" | tr ':' '\n'
#exit

descriptions=( \
"2-sample rRNAs, bowtie2" \
"2-sample rRNAs, STAR" \
)

args=( \
"--target 16S.fa 23S.fa
--out ROC_test_results
--verbose
--random-primer-len 9
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus" \
\
"--target 16S.fa 23S.fa
--out ROC_test_results
--verbose
--random-primer-len 9
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus
--star-aligner" \
\
)

total_count=0
pass_count=0

echo '[==========] ROC tests'
echo '[----------]'


for i in "${!descriptions[@]}"; do
#for i in 1; do

    echo "[ RUN      ] ROC test ${i}: "${descriptions[${i}]}

    total_count=$((total_count+1))

    arg="${args[${i}]}"

    # add name param if not provided already
    if [[ ${arg} != *"--name"* ]]; then
        arg="${arg} --name ROC_${i}"
    fi

    SECONDS=0

    shapemapper \
    ${arg} \
    --overwrite \
    --min-depth 100 \
    2>&1 \
    >shapemapper_temp/out.txt

    echo "Run time: ${SECONDS}s"

    both_pass=true

    # shapemapper-2.1.4:
    # bowtie2:
    # small subunit AUC: 0.728
    # large subunit AUC: 0.709
    # STAR:
    # small subunit AUC: 0.729
    # large subunit AUC: 0.710
    
    # lowered thresholds slightly
    # previously were .7291, .7085

    map=ROC_test_results/ROC_${i}_16S.map
    python3 \
      ${BASE_DIR}/internals/bin/area_under_ROC_curve.py \
      --map "${map}" \
      --ct "16S.ct" \
      --min-auc 0.728 \
      --name "small subunit"
    if [ $? != 0 ]; then
        both_pass=false
    fi

    map=ROC_test_results/ROC_${i}_23S.map
    python3 \
      ${BASE_DIR}/internals/bin/area_under_ROC_curve.py \
      --map "${map}" \
      --ct "23S.ct" \
      --min-auc 0.707 \
      --name "large subunit" 
    if [ $? != 0 ]; then
        both_pass=false
    fi

    if $both_pass; then
        pass_count=$((pass_count+1))
        echo "[       OK ] ROC test ${i}: "${descriptions[${i}]}
    else
        #cat shapemapper_temp/out.txt
        echo "[  FAILED  ] ROC test ${i}: "${descriptions[${i}]}
    fi


done

echo '[----------]'
echo '[==========]'
echo "[  PASSED  ] ${pass_count} tests for end-to-end pipeline success."

if [[ $pass_count != $total_count ]]; then
    echo "[  FAILED  ] $((total_count-pass_count)) tests."
fi

rm -rf core.*
if [[ $pass_count == $total_count ]]; then
    # remove intermediate files and logs if all runs succeeded
    rm -rf *shapemapper_log.txt    
    rm -rf shapemapper_temp
    rm -rf ROC_test_results
    :
fi
