#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"

SHAPEMAPPER_USE_SYSTEM_LIBS=false
SHAPEMAPPER_STATIC_BOOST=true
SHAPEMAPPER_DEBUG_SYMBOLS=false

usage="Build shapemapper binaries.
By default uses local libs and tools and statically links to boost libs.

Options:
 --use-system-libs     Use system libs and tools and dynamically link to boost libs. Default=false
 --static-boost        Force static linking to boost libs
 --debug-symbols       Build with debugging symbols. Default=false
"


force_static_boost=false

while [ "$1" != "" ]; do
    case "$1" in
        -h|--help)
            echo "${usage}"
            exit
            ;;
        --use-system-libs)
            SHAPEMAPPER_USE_SYSTEM_LIBS=true
            SHAPEMAPPER_STATIC_BOOST=false
            ;;
        --static-boost)
            force_static_boost=true
            ;;
        --debug-symbols)
            SHAPEMAPPER_DEBUG_SYMBOLS=true
            ;;
        *)
            echo "Unrecognized option \"${1}\"";
            echo "${usage}";
            exit 1;;
    esac
    shift
done

if ${force_static_boost}; then
    SHAPEMAPPER_STATIC_BOOST=true
fi

echo "SHAPEMAPPER_USE_SYSTEM_LIBS = ${SHAPEMAPPER_USE_SYSTEM_LIBS}"
# FIXME: not sure this flag is used anymore. respect or remove
echo "SHAPEMAPPER_STATIC_BOOST = ${SHAPEMAPPER_STATIC_BOOST}"
# FIXME: handle in cmakelists
echo "SHAPEMAPPER_DEBUG_SYMBOLS = ${SHAPEMAPPER_DEBUG_SYMBOLS}"

if [ "${SHAPEMAPPER_USE_SYSTEM_LIBS}" = false ]; then
    echo -e "\n\nWill attempt to build using local libs and tools. (Use \"--use-system-libs\" to disable).\n\n"
    source ${BASE_DIR}/internals/install/curl_paths.sh
    source ${BASE_DIR}/internals/install/dep_paths.sh
fi

cmake_args=""
if [ "${SHAPEMAPPER_DEBUG_SYMBOLS}" = true ]; then
    cmake_args="-DCMAKE_BUILD_TYPE=RelWithDebInfo"
fi

cd ${BASE_DIR}/internals

echo "Building shapemapper executables"
#rm -rf build
mkdir -p build
cd build
cmake ../cpp-src ${cmake_args}
make all
cp src/shapemapper_* ../bin/

shopt -s nullglob
test_binaries=( test/test_* )
# suppress file not found errors (in case test binaries were not built)
if (( ${#test_binaries[@]} )); then
    cp test/test_* ../bin/
fi

# FIXME: make sure miniconda python is installed
#        and in PATH before running this
#echo "Compiling RING-MaP Cython code"
#cd ${BASE_DIR}/bin/single_molecule
#python3.5 setup.py build_ext --inplace
