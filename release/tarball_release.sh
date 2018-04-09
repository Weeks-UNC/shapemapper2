#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)

THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"

source ${BASE_DIR}/install/dep_paths.sh

minimal=0 
activedev=0
while :; do
case $1 in
    --minimal)
        minimal=1
        ;;
    --include-active-dev)
        activedev=1
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

cd ${BASE_DIR}/..

DIRNAME=$(basename ${BASE_DIR})

VERSION="$(<${BASE_DIR}/release/version.txt)"

EXCLUDES="${DIRNAME}/.git \
${DIRNAME}/*.idea \
${DIRNAME}/build_deps \
${DIRNAME}/build \
${DIRNAME}/release/*.tar.gz \
${DIRNAME}/thirdparty/curl-src \
${DIRNAME}/thirdparty/pv* \
${DIRNAME}/build_thirdparty.log \
${DIRNAME}/*temp \
${DIRNAME}/*tmp \
${DIRNAME}/*results \
${DIRNAME}/*out"

if [ $activedev == 1 ]; then
    VERSION+="-active-dev"
else
    EXCLUDES+=" ${DIRNAME}/bin/single_molecule"
fi

if [ $minimal == 1 ]; then
    EXCLUDES+=" ${DIRNAME}/thirdparty"
    VERSION+="-minimal"
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

mv ${tarball_name} ${BASE_DIR}/release/
