#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

# Run end-to-end success tests for several common pipelines
# TODO: test for presence of non-empty output files
# FIXME: make it easier to terminate testing script with CTRL-C
# FIXME: still have occasional zombie procs if script terminated early with CTRL-C

CURRENT_DIR="$(pwd)"

# Find the parent folder of this script,
# resolving (possibly nested) symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$BASE_DIR/$SOURCE"
done
BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && cd .. && pwd )"



export PATH=${BASE_DIR}/bin:${BASE_DIR}:${PATH}

source ${BASE_DIR}/install/thirdparty_paths.sh

cd ${BASE_DIR}/test/data

mkdir -p shapemapper_temp

# TODO: rework how these arrays are laid out so test name, index, and args are nearby

descriptions=( \
"1 sample parallel SHAPE" \
"1 sample parallel SHAPE, sequence correction, 2 processors" \
"1 sample parallel SHAPE, STAR aligner, sequence correction, 2 processors" \
"1 sample parallel SHAPE, multiple folders" \
"1 sample parallel SHAPE, structured output" \
"1 sample parallel SHAPE explicit filenames" \
"1-sample parallel sequence correction" \
"3-sample parallel sequence correction and SHAPE" \
"2-sample serial sequence correction and SHAPE" \
"1-sample unpaired parallel sequence correction and SHAPE, render flowchart" \
"3-sample parallel SHAPE, output intermediate SAM files, parsed mutations, and mutation counts" \
"2-sample, 2-target parallel SHAPE, verbose" \
"2-sample, 2-target parallel SHAPE, --target-raw" \
"2-sample, 2-target parallel sequence correction and SHAPE" \
"2-sample, 2-target SHAPE, debug mutation output" \
"2-sample, 1-target SHAPE, debug mutation output" \
"2-sample, 1-target SHAPE, STAR aligner, debug mutation output" \
"2-sample, 2-target SHAPE, STAR aligner, debug mutation output" \
"1 sample, SHAPE and sequence correction (insert from reference)" \
"1 sample, SHAPE and sequence correction (deletion from reference)" \
)

args=( \
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPplus" \
\
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPplus
--correct-seq
--folder TPPminus
--nproc 2" \
\
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPplus
--correct-seq
--folder TPPminus
--nproc 2
--star-aligner" \
\
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPdenat TPPminus" \
\
"--target TPP.fa
--out end-to-end_test_results
--structured-output
--modified
--folder TPPplus" \
\
"--target TPP.fa
--out end-to-end_test_results
--modified
--R1 TPPplus/TPPplus_R1_aa.fastq.gz
--R2 TPPplus/TPPplus_R2_aa.fastq.gz" \
\
"--target TPP.fa
--out end-to-end_test_results
--correct-seq
--folder TPPminus" \
\
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPplus
--untreated
--folder TPPminus
--denatured
--folder TPPdenat
--correct-seq
--folder TPPminus" \
\
"--target TPP.fa
--out end-to-end_test_results
--serial
--modified
--folder TPPplus
--untreated
--folder TPPminus
--correct-seq
--folder TPPminus" \
\
"--target TPP.fa
--out end-to-end_test_results
--render-flowchart
--correct-seq
--unpaired-folder TPPminus
--modified
--unpaired-folder TPPplus" \
\
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPplus
--untreated
--folder TPPminus
--denatured
--folder TPPdenat
--output-aligned
--output-parsed
--output-counted" \
\
"--target 16S_trunc.fa 23S_trunc.fa
--out end-to-end_test_results
--verbose
--random-primer-len 9
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus" \
\
"--target 16S_trunc.fa
--target-raw >23S\n"\
"AGATAGCTGGTTCTCCCCGAAAGCTATTTAGGTAGCGCCTC\n"\
"GTGAATTCATCTCCGGGGGTAGAGCACTGTTTCGGCAAGGGGGTCATCCCGACTTACCAA\n"\
"CCCGATGCAAACTGCGAATACCGGAGAATGTTATCACGGGAGACACACGGCGGGTGCTAA\n"\
"CGTCCGTCGTGAAGAGGGAAACAACCCAGACCGCCAGCTA "\
"--out end-to-end_test_results
--verbose
--random-primer-len 9
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus" \
\
"--target 16S_trunc.fa 23S_trunc.fa
--out end-to-end_test_results
--random-primer-len 9
--correct-seq
--folder ribosome_minus
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus" \
\
"--target 16S_trunc.fa 23S_trunc.fa
--out end-to-end_test_results
--render-mutations
--random-primer-len 9
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus" \
\
"--target 16S_trunc.fa
--out end-to-end_test_results
--output-counted-mutations
--render-mutations
--random-primer-len 9
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus" \
\
"--target 16S_trunc.fa
--out end-to-end_test_results
--star-aligner
--output-counted-mutations
--render-mutations
--random-primer-len 9
--modified
--unpaired-folder ribosome_only_mapping_plus
--untreated
--unpaired-folder ribosome_only_mapping_minus" \
\
"--target 16S_trunc.fa 23S_trunc.fa
--out end-to-end_test_results
--output-counted-mutations
--star-aligner
--render-mutations
--random-primer-len 9
--modified
--unpaired-folder ribosome_only_mapping_plus
--untreated
--unpaired-folder ribosome_only_mapping_minus" \
\
"--target TPP_with_gap.fa
--out end-to-end_test_results
--correct-seq
--folder TPPminus
--modified
--folder TPPplus" \
\
"--target TPP_with_insert.fa
--out end-to-end_test_results
--correct-seq
--folder TPPminus
--modified
--folder TPPplus" \
\
)

# TODO: add this test back in once bowtie2 wrapper quotes -x index argument filename (PR submitted)
#"1 sample parallel SHAPE, spaces in name" \
#"--target TPP.fa
#--name 'end-to-end with spaces in name'
#--out end-to-end_test_results
#--modified
#--folder TPPplus" \
#\

total_count=0
pass_count=0

echo '[==========] End-to-end pipeline success tests'
echo '[----------]'

for i in "${!descriptions[@]}"; do
#for i in 20; do

    echo "[ RUN      ] End-to-end test ${i}: "${descriptions[${i}]}

    total_count=$((total_count+1))

    arg="${args[${i}]}"

    # add name param if not provided already
    if [[ ${arg} != *"--name"* ]]; then
        arg="${arg} --name end-to-end_${i}"
    fi

    SECONDS=0

    shapemapper \
    ${arg} \
    --overwrite \
    --min-depth 10 \
    2>&1 \
    >shapemapper_temp/out.txt

    echo "Run time: ${SECONDS}s"

    # just look at last line
    out=$(tail -n 1 shapemapper_temp/out.txt)
    out4=$(tail -n 3 shapemapper_temp/out.txt | head -n 1)

    # check if output indicates completed run (regardless of possible 
    # quality control check failures)
    if [[ "$out" == "ShapeMapper run successfully completed"* || \
          "$out4" == "ShapeMapper run completed"* ]]; then
        pass_count=$((pass_count+1))
        echo "[       OK ] End-to-end test ${i}: "${descriptions[${i}]}
    else
        cat shapemapper_temp/out.txt
        echo "[  FAILED  ] End-to-end test ${i}: "${descriptions[${i}]}
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
    rm -rf end-to-end_test_results
    :
fi
