#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
BUILD_DEPS=${BASE_DIR}/build_deps

cd ${BASE_DIR}

# Check for curl, build if necessary
#export PATH=${BUILD_DEPS}/curl-7.47.1/bin:$PATH
# Check for curl, unzip and build if not present
#v=$(curl --version)
#set ${v}
#char=${2:0:1}
#FOUND_CURL=false
#if [ "${char}" -eq "${char}" ] 2>/dev/null; then
#  # is a numeric char
#  echo "Curl already present, version is ${2}"
#  FOUND_CURL=true
#else
#  # is not a numeric char
#  echo "curl not found, will attempt to build"
#fi


# need dev headers and lib for cmake build, so default system curl
# usually isn't enough. Just build it.

#if [ "${FOUND_CURL}" = false ]; then
if [ ! -f "${BUILD_DEPS}/curl-7.47.1/bin/curl" ]; then
    echo "Building curl..."
    mkdir -p thirdparty
    cd thirdparty
    rm -rf curl-7.47.1
    rm -rf curl-src
    tar -xf ${BASE_DIR}/install/curl-7.47.1.tar.gz
    mv curl-7.47.1 curl-src
    cd curl-src
    make prefix=${BUILD_DEPS}/curl-7.47.1
    make prefix=${BUILD_DEPS}/curl-7.47.1 install
fi


