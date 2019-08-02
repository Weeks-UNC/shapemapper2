#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)

THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"


source ${BASE_DIR}/internals/paths/dep_paths.sh

sourceonly=0
while :; do
case $1 in
    --source-only)
        sourceonly=1
        ;;
    -?*)
        printf 'ERROR: Unknown argument: %s\n' "$1" >&2
        exit 1
        ;;
    *)
        break
esac
shift
done


TARBALL_DIR=$(readlink -f ${BASE_DIR}/..)
DIRNAME=$(basename ${BASE_DIR})

cd ${BASE_DIR}/..

VERSION="$(<${BASE_DIR}/internals/release/version.txt)"

EXCLUDES="${DIRNAME}/.git \
${DIRNAME}/*.idea \
${DIRNAME}/internals/build_deps \
${DIRNAME}/build \
${DIRNAME}/internals/build \
${DIRNAME}/internals/release/*.tar.gz \
${DIRNAME}/internals/thirdparty/*.tar.gz \
${DIRNAME}/internals/thirdparty/*.zip \
${DIRNAME}/internals/thirdparty/curl-src \
${DIRNAME}/internals/thirdparty/pv* \
${DIRNAME}/internals/build_thirdparty.log \
${DIRNAME}/internals/python/run_debug.py \
${DIRNAME}/internals/*temp \
${DIRNAME}/internals/*tmp \
${DIRNAME}/internals/*results \
${DIRNAME}/*out"


if [ $sourceonly == 1 ]; then
    VERSION+="-source-only"
    EXCLUDES+=" ${DIRNAME}/internals/thirdparty \
${DIRNAME}/internals/bin/shapemapper_read_trimmer \
${DIRNAME}/internals/bin/shapemapper_mutation_counter \
${DIRNAME}/internals/bin/shapemapper_mutation_parser \
${DIRNAME}/internals/bin/test_histogram \
${DIRNAME}/internals/bin/test_mutation_counter \
${DIRNAME}/internals/bin/test_mutation_parser \
${DIRNAME}/internals/bin/test_read_trimmer"
fi

tarball_name="shapemapper-${VERSION}.tar.gz"
tarball_path="${TARBALL_DIR}"

s=""
for x in ${EXCLUDES}; do
    s=${s}" --exclude=${x}"
done

s=${s}" --exclude='*.out'"

# rename top-level folder to indicate release
# (this results in redundant path separators)
s="--transform=s,${DIRNAME},shapemapper-${VERSION}/, $s"

cmd="tar -cvpzf ${tarball_name} --exclude-backups --exclude-vcs --show-transformed-names $s ${DIRNAME}"
echo "$cmd"

$cmd

mv ${tarball_name} ${BASE_DIR}/internals/release/
