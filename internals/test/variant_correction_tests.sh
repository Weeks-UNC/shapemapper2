#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

# Run sequence variant correction tests

# FIXME: make it easier to terminate testing script with CTRL-C
# TODO: add pipeline latency options to CLI for faster test execution
#       (already have param in pipeline and module run() methods)

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

descriptions=("unchanged" \
"ambig del G near 5-prime end" \
"ambig del G near 5-prime at different location" \
"ambig del G at high-background position" \
"unambig del C near 5-prime" \
"ambig del T near 5-prime" \
"mismatch G->T near 5-prime" \
"double mismatch GC->AT near 5-prime" \
"insert CCA towards 5-prime end" \
"GT->C toward 5-prime" \
"GT->C toward 5-prime, del T nearby, GT->CA toward 3-prime" \
"mismatch G->T next to del A" \
"mismatch T->C next to insert A" \
)

vars=( \
"GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGGGCCAAGGACTCGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGTGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGGATCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGGGCCAAGGACTCGGGGCCATGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGGGCCAAGGACTCGGGCGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGGGCCAAGGACTCGGGCGCCCTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGCATCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGGGCCAAGTCTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
"GGCCTTCGGGCCAAGGACTCGGGGCAGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC" \
)

total_count=0
pass_count=0

echo '[==========] Sequence variant correction'
echo '[----------]'

for i in "${!descriptions[@]}"; do
    v=${vars[${i}]}
    echo "[ RUN      ] Sequence variant ${i}: "${descriptions[${i}]}

    total_count=$((total_count+1))

    echo -e ">TPP\n"${v} > "shapemapper_temp/variant_test.fa"

    shapemapper \
    --target shapemapper_temp/variant_test.fa \
    --out variant_test_results \
    --name variant_${i} \
    --correct-seq \
    --min-seq-depth 40 \
    --overwrite \
    --folder TPPminus >/dev/null

    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

    out=variant_test_results/variant_${i}*.fa

    python3 ${BASE_DIR}/internals/bin/compare_fasta.py "TPP.fa" $out \
    2> shapemapper_temp/stderr.txt

    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

    stderr=$(<shapemapper_temp/stderr.txt)

    if [[ "$stderr" != ERROR* ]]; then
        pass_count=$((pass_count+1))
        echo "[       OK ] Sequence variant ${i}: "${descriptions[${i}]}
    else
        echo $stderr
        echo "[  FAILED  ] Sequence variant ${i}: "${descriptions[${i}]}
    fi

done

echo '[----------]'
echo '[==========]'
echo "[  PASSED  ] ${pass_count} tests for sequence variant correction."

if [[ $pass_count != $total_count ]]; then
    echo "[  FAILED  ] $((total_count-pass_count)) tests."
fi

rm -rf core.*
if [[ $pass_count == $total_count ]]; then
    # remove intermediate files and logs if all runs succeeded
    rm -rf *shapemapper_log.txt
    rm -rf shapemapper_temp
    rm -rf variant_test_results
    :
fi
