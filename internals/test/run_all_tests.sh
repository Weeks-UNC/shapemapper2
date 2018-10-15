#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

# Run all tests:
#  - c++ unit tests
#  - end-to-end pipeline tests for success
#  - end-to-end pipeline tests for specific component failure detection
#  - end-to-end sequence variant correction tests

# TODO: suppress "tput: No value for $TERM" warnings


# Find the parent folder of this script,
# resolving (possibly nested) symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    TEST_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$TEST_DIR/$SOURCE"
done
TEST_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

tests=( \
"${TEST_DIR}/cpp_unit_tests.sh" \
"${TEST_DIR}/end-to-end_tests.sh" \
"${TEST_DIR}/component_failure_tests.sh" \
"${TEST_DIR}/variant_correction_tests.sh" \
"${TEST_DIR}/ROC_tests.sh" \
)

names=( \
"c++ unit" \
"end-to-end success" \
"module failure detection" \
"sequence variant correction" \
"area under ROC curve" \
)

overall_total_count=0
overall_fail_count=0
total_counts=( 0 0 0 0 0 )
fail_counts=( 0 0 0 0 0 )

exec 5>&1

for i in "${!tests[@]}"; do
    test=${tests[${i}]}
    name=${names[${i}]}

    echo -e "Running ${name} tests . . .\n"

    # want to display script output while running,
    # and capture output in variable, and give error if failed
    # FIXME: Component failure test script output doesn't flush except on error
    out=$($test 2>&1 | tee >(cat - >&5); echo ${PIPESTATUS[0]})
    # - Tricky to also get return code cuz subshell + pipes

    # get return code
    rc=$(echo "${out}" | tail -n 1)
    # strip return code from string
    out=$(echo "${out}" | head -n -1)
    if [[ $rc != 0 ]]; then
        echo "ERROR: a test script itself failed."
        exit $rc
    fi

    # tests that "RUN" but don't "OK" have failed - identify these
    # by matching RUN lines to OK lines
    fail_count=0
    total_count=0
    run_lines=$(echo "${out}" | grep '\[ RUN      ]')
    ok_lines=$(echo "${out}" | grep '\[       OK ]')
    failed_lines=""
    while read -r run_line; do
        total_count=$((total_count+1))
        stripped_run=$(echo "$run_line" | cut -c 14-999 | sed 's/^ *//g' | sed 's/ *$//g')
        found_matching_ok=false
        while read -r ok_line; do
            stripped_ok=$(echo "$ok_line" | cut -c 14-999 | sed 's/^ *//g')
            if [[ $stripped_ok == $stripped_run* ]]; then
                found_matching_ok=true
                break
            fi
        done <<< "$ok_lines"
        if [[ $found_matching_ok != true ]]; then
            #echo "FAILED test found: $stripped_run"
            fail_count=$((fail_count+1))
        fi
    done <<< "$run_lines"
    total_counts[${i}]=$total_count
    fail_counts[${i}]=$fail_count
    overall_total_count=$((overall_total_count+total_count))
    overall_fail_count=$((overall_fail_count+fail_count))
    # print line the width of terminal
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
    echo -e "\n\n"
done


if [[ $overall_fail_count == 0 ]]; then
    echo "All tests passed"
    echo "SUCCESS"
    exit 0
else
    # summarize tests
    for i in 0 1 2 3 4; do
        test=${tests[${i}]}
        name=${names[${i}]}
        total_count=${total_counts[${i}]}
        fail_count=${fail_counts[${i}]}
        echo -e "${fail_count} / ${total_count}\t${name} test(s) failed."
    done
    echo "----------"
    echo -e "${overall_fail_count} / ${overall_total_count}\ttotal test(s) failed."
    echo "FAILURE"
    exit 1
fi
