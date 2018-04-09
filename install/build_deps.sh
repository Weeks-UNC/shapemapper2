#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
BUILD_DEPS=${BASE_DIR}/build_deps

# copy all output to logfile
#exec > >(tee "$BASE_DIR/build_deps.log") 2>&1

source ${BASE_DIR}/install/curl_paths.sh

cd ${BASE_DIR}

mkdir -p ${BUILD_DEPS}

# install curl
${BASE_DIR}/install/build_curl.sh

cd ${BUILD_DEPS}

# gcc, g++ (older versions don't support C11 features)
if [ ! -f "gcc/bin/gcc" ]; then
    #p='http://ftp.gnu.org/gnu/gsrc/gsrc-2014.10.11.tar.gz' # This package could be easier to work with
    p='http://ftp.gnu.org/gnu/gcc/gcc-5.3.0/gcc-5.3.0.tar.gz'
    echo "Downloading gcc source"
    curl -L $p -o gcc-5.3.0.tar.gz
    rm -rf gcc-5.3.0
    rm -rf gcc-src
    echo "Extracting gcc source archive"
    tar -xf gcc-5.3.0.tar.gz    
    mv gcc-5.3.0 gcc-src
    echo "Configuring gcc build"
    cd gcc-src
    contrib/download_prerequisites
    cd ${BUILD_DEPS}
    rm -rf gcc-build
    rm -rf ${BUILD_DEPS}/gcc
    mkdir gcc-build
    cd gcc-build
    echo "LIBRARY_PATH: "${LIBRARY_PATH}
    ../gcc-src/configure --prefix=${BUILD_DEPS}/gcc --disable-multilib --enable-languages=c,c++
    echo "Building gcc"
    cd ${BUILD_DEPS}/gcc-build
    make
    make install
    # symlink using generic name so cmake will find the correct compiler
    cd ${BUILD_DEPS}/gcc/bin
    ln -s gcc cc
fi

cd ${BASE_DIR}

source ${BASE_DIR}/install/dep_paths.sh

cd ${BUILD_DEPS}

#echo "about to build git"
#echo "LIBRARY_PATH: ${LIBRARY_PATH}"
#echo "CPATH: ${CPATH}"
#exit

# bare-bones git (needed for cmake automagic download features)
if [ ! -f "git-2.7.2/bin/git" ]; then
    p='https://github.com/git/git/tarball/v2.7.2'
    echo "Downloading git source"
    curl -L $p -o git-2.7.2.tar.gz
    echo "Extracting git source archive"
    rm -rf git-2.7.2
    rm -rf git-src
    tar -xf git-2.7.2.tar.gz
    mv git-git-e4980cd git-src
    echo "Building git"
    cd git-src
    make prefix=${BUILD_DEPS}/git-2.7.2 \
     NO_PERL=YesPlease \
     NO_TCLTK=YesPlease \
     all 
    make prefix=${BUILD_DEPS}/git-2.7.2 \
     NO_PERL=YesPlease \
     NO_TCLTK=YesPlease \
     install
    cd ${BUILD_DEPS}
fi

cd ${BUILD_DEPS}

# cmake
if [ ! -f "cmake-3.4.3-Linux-x86_64/bin/cmake" ]; then
    p='https://cmake.org/files/v3.4/cmake-3.4.3-Linux-x86_64.tar.gz'
    echo "Downloading cmake"
    curl --insecure -L $p -o cmake-3.4.3.tar.gz # certificate problems for cmake.org, so use --insecure
    echo "Extracting cmake archive"
    tar -xf cmake-3.4.3.tar.gz
fi

cd ${BUILD_DEPS}

# newer version of tar (for releases)
# - old versions don't support --exclude-backups, other nice options
# - Currently having problems building on killdevil with gcc 5.3
if [ ! -f "tar/bin/tar"  ]; then
    p='http://ftp.gnu.org/gnu/tar/tar-1.28.tar.gz'
    echo "Downloading tar source"
    curl -L $p -o tar-1.28.tar.gz
    rm -rf tar-1.28
    rm -rf tar-src
    echo "Extracting tar source archive"
    tar -xf tar-1.28.tar.gz
    mv tar-1.28 tar-src
    rm -rf tar-build
    rm -rf ${BUILD_DEPS}/tar
    mkdir tar-build
    cd tar-build
    export CFLAGS='-gdwarf-2 -gstrict-dwarf'
    ../tar-src/configure --prefix=${BUILD_DEPS}/tar
    echo "Building tar"
    cd ${BUILD_DEPS}/tar-build
    make
    make install
fi

cd ${BUILD_DEPS}

# zlib
#if [ ! -f zlib/lib/libz.a ]; then
#    p='http://zlib.net/zlib-1.2.8.tar.gz'
#    echo "Downloading zlib"
#    curl -L ${p} -o zlib-1.2.8.tar.gz
#    echo "Extracting zlib archive"
#    tar -xf zlib-1.2.8.tar.gz
#    cd zlib-1.2.8
#    ./configure --static --prefix=${BUILD_DEPS}/zlib
#    make
#    make install
#fi
#
#cd ${BUILD_DEPS}

# boost
boost_libs="filesystem program_options iostreams system"
found_all_libs=true
for lib in ${boost_libs}; do
    if [ ! -f boost/lib/libboost_${lib}.a ]; then
        found_all_libs=false
        break
    fi
done

if [ ${found_all_libs} = true ]; then
    echo "All static boost libs already built."
else
    rm -rf boost
    rm -rf boost-build
    rm -rf boost_1_60_0
    rm -rf boost_1_60_0-src

    # full boost source
    #p='http://iweb.dl.sourceforge.net/project/boost/boost/1.60.0/boost_1_60_0.tar.gz'
    p='http://downloads.sourceforge.net/project/boost/boost/1.60.0/boost_1_60_0.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fboost%2Ffiles%2Fboost%2F1.60.0%2F&ts=1458926317&use_mirror=iweb'

    echo "Downloading Boost source"
    curl -L ${p} -o boost_1_60_0.tar.gz

    echo "Extracting Boost source archive"
    tar -xf boost_1_60_0.tar.gz
    mv boost_1_60_0 boost_1_60_0-src

    pre=${BUILD_DEPS}/boost
    build=${BUILD_DEPS}/boost-build
    cd boost_1_60_0-src

    echo "Building Boost build system"
    ./bootstrap.sh

    for lib in ${boost_libs}; do
        if [ ! -f ${pre}/lib/libboost_${lib}.a ]; then
            echo "Building Boost ${lib} library"
            ./b2 install --prefix=${pre} --build-dir=${build} --with-${lib} link=static
        else
            echo "Boost $lib library already built, skipping."
        fi
    done
fi

## Doxygen and Python3.5 + Cython + numpy (install through conda)
## - need python to compile RING-MaP cython code

# install miniconda if not present
${BASE_DIR}/install/get_miniconda.sh # This script depends on curl
source ${BASE_DIR}/install/dep_paths.sh

unset PYTHONPATH
CONDA_PATH=${BASE_DIR}/thirdparty/miniconda
export PATH=${CONDA_PATH}/bin:$PATH

#echo -en "y\n" | conda install -c bioconda cython numpy
echo -en "y\n" | conda install -c omnia doxygen
# create a link to doxygen executable so we can add it to PATH
# later without loading everything else in the miniconda bin
ln -sf ${CONDA_PATH}/bin/doxygen ${BUILD_DEPS}/doxygen
# also create a link to python3.5
#ln -sf ${CONDA_PATH}/bin/python3.5 ${BUILD_DEPS}/python3.5
