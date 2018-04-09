#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
BUILD_DEPS=${BASE_DIR}/build_deps

# curl
PATH=${BUILD_DEPS}/curl-7.47.1/bin:${PATH}
CPATH=${BUILD_DEPS}/curl-7.47.1/include:${CPATH}
LIBRARY_PATH=${BUILD_DEPS}/curl-7.47.1/lib:${LIBRARY_PATH}
LD_LIBRARY_PATH=${BUILD_DEPS}/curl-7.47.1/lib:${LD_LIBRARY_PATH}

# Strip trailing colon if present and export
# (otherwise a blank initial LIBRARY_PATH will mess up gcc config)
for v in PATH CPATH LIBRARY_PATH LD_LIBRARY_PATH; do
    [ "${!v: -1}" == ":" ] && tmp=${!v%?} && eval $v=\${tmp}
    export $v
done
