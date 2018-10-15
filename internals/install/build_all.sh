#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"

if [[ "$@" != "--use-system-libs" ]]; then
    # download and build local libs and tools if not using system libs
    ${BASE_DIR}/internals/install/build_deps.sh
fi

${BASE_DIR}/internals/install/build_binaries.sh "${@}"

