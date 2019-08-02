#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"
BUILD_DEPS=${BASE_DIR}/internals/build_deps

# gcc/g++
PATH=${BUILD_DEPS}/gcc/bin:${PATH}
LIBRARY_PATH=${BUILD_DEPS}/gcc/lib64:${LIBRARY_PATH}
LD_LIBRARY_PATH=${BUILD_DEPS}/gcc/lib64:${LD_LIBRARY_PATH}

# cmake
PATH=${BUILD_DEPS}/cmake-3.4.3-Linux-x86_64/bin:${PATH}

# git
PATH=${BUILD_DEPS}/git-2.7.2/bin:${PATH}

# tar
PATH=${BUILD_DEPS}/tar/bin:${PATH}

# zlib
#CPATH=${BUILD_DEPS}/zlib/include:${CPATH}
#LIBRARY_PATH=${BUILD_DEPS}/zlib/lib:${LIBRARY_PATH}

# boost
BOOST_ROOT=${BUILD_DEPS}/boost
BOOST_INCLUDEDIR=${BUILD_DEPS}/boost/include
BOOST_LIBRARYDIR=${BUILD_DEPS}/boost/lib
CMAKE_ARGS="-DBoost_USE_STATIC_LIBS=ON "
CMAKE_ARGS="-DBOOST_ROOT=${BOOST_ROOT} ${CMAKE_ARGS}"
CMAKE_ARGS="-DBOOST_INCLUDEDIR=${BOOST_INCLUDEDIR} ${CMAKE_ARGS}"
CMAKE_ARGS="-DBOOST_LIBRARYDIR=${BOOST_LIBRARYDIR} ${CMAKE_ARGS}"

# Doxygen (there should be a symlink to thirdparty/miniconda/bin/doxygen)
# Python (there should be a symlink to thirdparty/miniconda/bin/python3.5)
PATH=${BUILD_DEPS}:${PATH}


# Strip trailing colon if present and export
# (otherwise a blank initial LIBRARY_PATH will mess up gcc config)
for v in PATH CPATH LIBRARY_PATH LD_LIBRARY_PATH; do
    [ "${!v: -1}" == ":" ] && tmp=${!v%?} && eval $v=\${tmp}
    export $v
done
