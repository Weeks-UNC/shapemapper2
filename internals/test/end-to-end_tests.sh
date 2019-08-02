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
BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && cd ../.. && pwd )"


source ${BASE_DIR}/internals/paths/bin_paths.sh

cd ${BASE_DIR}/internals/test/data

mkdir -p shapemapper_temp


names_args_interleaved=( \
"1 sample parallel SHAPE" \
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPplus" \
\
"1 sample parallel SHAPE, unpaired input" \
"--target TPP.fa
--out end-to-end_test_results
--modified
--unpaired-folder TPPplus" \
\
"1 sample parallel SHAPE, sequence correction, 2 processors" \
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPplus
--correct-seq
--folder TPPminus
--nproc 2" \
\
"1 sample parallel SHAPE, STAR aligner, sequence correction, 2 processors" \
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPplus
--correct-seq
--folder TPPminus
--nproc 2
--star-aligner" \
\
"1 sample parallel SHAPE, STAR aligner, explicit genomeSAindexNbase" \
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPplus
--genomeSAindexNbase 3
--star-aligner" \
\
"1 sample parallel SHAPE, multiple folders" \
"--target TPP.fa
--out end-to-end_test_results
--modified
--folder TPPdenat TPPminus" \
\
"1 sample parallel SHAPE, multiple folders, unpaired input" \
"--target TPP.fa
--out end-to-end_test_results
--modified
--unpaired-folder TPPdenat TPPminus" \
\
"1 sample parallel SHAPE, multiple folders, unpaired input, STAR aligner" \
"--target TPP.fa
--out end-to-end_test_results
--modified
--star-aligner
--unpaired-folder TPPdenat TPPminus" \
\
"1 sample parallel SHAPE, structured output" \
"--target TPP.fa
--out end-to-end_test_results
--structured-output
--modified
--folder TPPplus" \
\
"1 sample parallel SHAPE explicit filenames" \
"--target TPP.fa
--out end-to-end_test_results
--modified
--R1 TPPplus/TPPplus_R1_aa.fastq.gz
--R2 TPPplus/TPPplus_R2_aa.fastq.gz" \
\
"1-sample parallel sequence correction" \
"--target TPP.fa
--out end-to-end_test_results
--correct-seq
--folder TPPminus" \
\
"3-sample parallel sequence correction and SHAPE" \
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
"2-sample serial sequence correction and SHAPE" \
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
"1-sample unpaired parallel sequence correction and SHAPE, render flowchart" \
"--target TPP.fa
--out end-to-end_test_results
--render-flowchart
--correct-seq
--unpaired-folder TPPminus
--modified
--unpaired-folder TPPplus" \
\
"3-sample parallel SHAPE, output intermediate SAM files, parsed mutations, and mutation counts" \
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
"2-sample, 2-target parallel SHAPE, verbose" \
"--target 16S_trunc.fa 23S_trunc.fa
--out end-to-end_test_results
--verbose
--random-primer-len 9
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus" \
\
"2-sample, 2-target parallel SHAPE, --target-raw" \
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
"2-sample, 2-target parallel sequence correction and SHAPE" \
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
"2-sample, 2-target SHAPE, debug mutation output" \
"--target 16S_trunc.fa 23S_trunc.fa
--out end-to-end_test_results
--render-mutations
--random-primer-len 9
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus" \
\
"2-sample, 1-target SHAPE, debug mutation output" \
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
"2-sample, 1-target SHAPE, STAR aligner, debug mutation output" \
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
"2-sample, 2-target SHAPE, STAR aligner, debug mutation output" \
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
"1 sample, SHAPE and sequence correction (insert from reference)" \
"--target TPP_with_gap.fa
--out end-to-end_test_results
--correct-seq
--folder TPPminus
--modified
--folder TPPplus" \
\
"1 sample, SHAPE and sequence correction (deletion from reference)" \
"--target TPP_with_insert.fa
--out end-to-end_test_results
--correct-seq
--folder TPPminus
--modified
--folder TPPplus" \
\
"amplicon primer pair mapping filter, primers located in sequence" \
"--target primer_pairs/RMRP/Rmrp.fa \
--out end-to-end_test_results
--overwrite
--min-depth 1000
--amplicon
--modified --folder primer_pairs/RMRP/plus" \
\
"amplicon primer pair mapping filter, primers in file" \
"--target primer_pairs/RMRP/Rmrp.fa \
--out end-to-end_test_results
--overwrite
--min-depth 1000
--amplicon
--primers primer_pairs/RMRP/primers_no_header.txt
--modified --folder primer_pairs/RMRP/plus" \
\
"amplicon primer pair mapping filter, STAR aligner" \
"--target primer_pairs/RMRP/Rmrp.fa \
--out end-to-end_test_results
--overwrite
--min-depth 1000
--amplicon
--star-aligner
--modified --folder primer_pairs/RMRP/plus" \
\
"STAR rerun after segfault" \
"--target TPP_masked_primers.fa
--out example_results
--log example_results/log.txt
--name TPP
--min-depth 1000
--modified --folder TPPplus
--amplicon
--star-aligner
--genomeSAindexNbase 5
--rerun-on-star-segfault
--rerun-genomeSAindexNbase 3" \
\
"2-sample, 2-target parallel SHAPE, bowtie2 effort params" \
"--target 16S_trunc.fa 23S_trunc.fa
--out end-to-end_test_results
--max-search-depth 12
--max-reseed 1
--random-primer-len 9
--modified
--folder ribosome_plus
--untreated
--folder ribosome_minus" \
\
)

# 29 total tests


# associative arrays not available in bash 3
# - abusing aliases to simulate

collapse() {
    echo "$1" | tr $'\n' ' ' | tr -s ' '
}

# map_set <name> <key> <value>
map_set() {
    alias "${1}$2"="$(collapse "$3")"
}

# map_get <name> <key>
map_get() {
    alias "${1}$2" | awk -F"'" '{ print $2; }'
}


# fill simulated array with key "row_col": value, 
# where 1st column is test name, 2nd column is arguments
len=$((${#names_args_interleaved[@]} / 2))
for ((i=0;i<len;i++)); do
    n=$((i*2))
    n1=$(($n+1))

    #echo "##################################"
    name="${names_args_interleaved[${n}]}"
    args="${names_args_interleaved[${n1}]}"
    #echo "${name}"
    #echo " - - - - - - - - "
    #echo "${args}"
    map_set a "${i}_0" "$name"
    map_set a "${i}_1" "$args"
done


# TODO: add this test back in once bowtie2 wrapper quotes -x index argument filename (PR submitted)
#"1 sample parallel SHAPE, spaces in name" \
#"--target TPP.fa
#--name 'end-to-end with spaces in name'
#--out end-to-end_test_results
#--modified
#--folder TPPplus" \
#

total_count=0
pass_count=0

echo '[==========] End-to-end pipeline success tests'
echo '[----------]'

for ((i=0;i<len;i++)); do
#for i in 0; do

    description=$(map_get a "${i}_0")
    arg=$(map_get a "${i}_1")


    echo "[ RUN      ] End-to-end test ${i}: "${description}

    total_count=$((total_count+1))

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

    out=$(tail -n 20 shapemapper_temp/out.txt)

    # strip warning message lines
    out=$(sed '/WARNING:/d' <<< "${out}")
    out=$(sed '/See quality control checks/d' <<< "${out}")

    # look at last line
    out=$(tail -n 1 <<< "${out}")

    # check if output indicates completed run (regardless of possible
    # quality control check failures)
    if [[ "$out" == "ShapeMapper run successfully completed"* || \
          "$out" == "ShapeMapper run completed"* ]]; then
        pass_count=$((pass_count+1))
        echo "[       OK ] End-to-end test ${i}: "${description}
    else
        cat shapemapper_temp/out.txt
        echo "[  FAILED  ] End-to-end test ${i}: "${description}
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
