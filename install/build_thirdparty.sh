#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2017 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"

cd ${BASE_DIR}

# copy all output to logfile
exec &> >(tee "$BASE_DIR/build_thirdparty.log")

# install curl and/or miniconda if not present
${BASE_DIR}/install/build_curl.sh
${BASE_DIR}/install/get_miniconda.sh # This script depends on curl

#source ${BASE_DIR}/install/curl_paths.sh # this is actually redundant, since already sourced in get_miniconda
source ${BASE_DIR}/install/dep_paths.sh

unset PYTHONPATH # not sure if this is working - conda still gives a warning
CONDA_PATH=${BASE_DIR}/thirdparty/miniconda
export PATH=${CONDA_PATH}/bin:$PATH

jdk_version=8.0.45-0
bbmap_version=37.52

# bowtie2 wrapper requires perl, which doesn't appear to be an explicit dependency in the bioconda package
conda install -yq \
-c bioconda \
"python=3.5.3" \
"perl-threaded=5.22.0" \
"bowtie2=2.2.7" \
"bbmap=${bbmap_version}" \
"java-jdk=${jdk_version}" \
"star=2.5.2a" \
"matplotlib=1.5.1"

# - need a version of graphviz built with old GLIBC for older platforms
# - Working graphviz=2.38.0 from anaconda channel uses GLIBC_2.3, GLIBC_2.3.4, and GLIBC_2.2.5
# This graphviz gives a warning, even though the linked file is apparently present
#   Warning: Could not load "<whatever>/thirdparty/miniconda/lib/graphviz/libgvplugin_pango.so.6" - file not found

conda install -yq \
-c anaconda \
"graphviz=2.38.0"

# graphviz from anaconda channel appears to pull in an incompatible libstdc++.so
# (breaking bowtie2), but graphviz still seems to work fine if the offending
# package is removed
conda remove -yq --force libstdcxx-ng

# pipeviewer                                                                                            
cd ${BASE_DIR}/thirdparty
p='http://www.ivarch.com/programs/sources/pv-1.6.0.tar.gz'
if [ ! -f 'pv-1.6.0.tar.gz' ]; then
    echo "Downloading Pipe Viewer utility"
    curl -L ${p} -o pv-1.6.0.tar.gz
fi
if [ ! -d pv-1.6.0 ]; then
    tar -xf pv-1.6.0.tar.gz
fi
if [ ! -d 'pipeviewer' ]; then
    cd pv-1.6.0
    ./configure --prefix=${BASE_DIR}/thirdparty/pipeviewer
    make
    make install
fi

# build bbmerge JNI components
export JAVA_HOME=${CONDA_PATH}/pkgs/java-jdk-${jdk_version}
cd ${CONDA_PATH}/opt/bbmap-${bbmap_version}/jni/
make -f makefile.linux


# check if environment appears to load correctly
# FIXME: move this to a more comprehensive test suite to run post-installation

set -x

source "${BASE_DIR}/install/thirdparty_paths.sh"

which bowtie2-build
bowtie2-build --version
which bowtie2
bowtie2 --version
which bbmerge.sh
bbmerge.sh --version
which pv
pv --version
which python3.5
python3.5 --version
which dot
dot -V



