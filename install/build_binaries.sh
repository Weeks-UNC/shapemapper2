#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"

SHAPEMAPPER_USE_SYSTEM_LIBS=false
SHAPEMAPPER_STATIC_BOOST=true

usage="Build shapemapper binaries.
By default uses local libs and tools and statically link to boost libs.

Options:
 --use-system-libs     Use system libs and tools and dynamically link to boost libs. Default=false
 --static-boost        Force static linking to boost libs
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
echo "SHAPEMAPPER_STATIC_BOOST = ${SHAPEMAPPER_STATIC_BOOST}"

if [ "${SHAPEMAPPER_USE_SYSTEM_LIBS}" = false ]; then
    echo -e "\n\nWill attempt to build using local libs and tools. (Use \"--use-system-libs\" to disable).\n\n"
    source ${BASE_DIR}/install/curl_paths.sh
    source ${BASE_DIR}/install/dep_paths.sh
fi


cd ${BASE_DIR}

echo "Building shapemapper executables"
#rm -rf build
mkdir -p build
cd build
cmake ../cpp-src ${CMAKE_ARGS}
make all
cp src/shapemapper_* ../bin/
cp test/test_* ../bin/

# FIXME: make sure miniconda python is installed
#        and in PATH before running this
#echo "Compiling RING-MaP Cython code"
#cd ${BASE_DIR}/bin/single_molecule
#python3.5 setup.py build_ext --inplace
