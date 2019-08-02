#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

set -e # exit on first error (if any)
set -x
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"

cd ${BASE_DIR}/internals

# copy all output to logfile
exec &> >(tee "$BASE_DIR/internals/build_thirdparty.log")

# install miniconda if not present
${BASE_DIR}/internals/thirdparty_helper/get_miniconda.sh

source ${BASE_DIR}/internals/paths/dep_paths.sh

unset PYTHONPATH # not sure if this is working - conda still gives a warning
CONDA_PATH=${BASE_DIR}/internals/thirdparty/miniconda
export PATH=${CONDA_PATH}/bin:$PATH

test_matplotlib () {
    python -c $'import matplotlib as mp\nmp.use("Agg")\nprint(mp.__version__)\nfrom matplotlib import pyplot as plt\nprint("successfully imported pyplot")'
}

jdk_version=8.0.45-0
bbmap_version=37.78 # 37.52
# recent bowtie2 release builds (not through conda) seem to be portable, so just use that

# note: need the --env option to prevent config option being written outside the thirdparty directory
# without --env option, conda config stuff is put in /<user>/.conda and /<user>/.condarc (user is "root" on buildbox)

# set channel priorities
conda config --env --add channels defaults
conda config --env --add channels conda-forge
conda config --env --add channels bioconda
#conda config --env --add channels anaconda # testing


conda config --env --set auto_update_conda false # updating conda itself sometimes requires newer GLIBC, so disable
conda config --env --set always_yes true

# bowtie2 wrapper requires perl

# need python=3.5.3 from defaults,
# bbmap, java, star, perl from bioconda,
# matplotlib from ?
# graphviz from anaconda

# matplotlib 2.1.0 brings in
#libstdcxx-ng:                  7.2.0 conda-forge
# which breaks conda bowtie2
# (probably related to collisions between libgcc and libstdcxx-ng)
# graphviz pulls in an incompatible libstdc++.so, but
# reinstalling libgcc after seems to work

# working perl: perl-threaded-5.22.0-10 from bioconda
# nonworking perl: perl-threaded-5.22.0-13 from bioconda requires GLIBC_2.11

# - Working graphviz=2.38.0 from anaconda channel uses GLIBC_2.3, GLIBC_2.3.4, and GLIBC_2.2.5
#   - gives a warning, even though the linked file is apparently present
#   Warning: Could not load "<whatever>/thirdparty/miniconda/lib/graphviz/libgvplugin_pango.so.6" - file not found


echo -e "\n### Installing main packages###\n"

# okay now getting same error with this command, so I'm pretty sure there's a server glitch somewhere

conda install -y \
-c defaults \
"python=3.5.3" \
"libgcc=5.2.0" \
"bbmap=${bbmap_version}" \
"java-jdk=${jdk_version}" \
"perl-threaded=5.22.0=10" \
"star=2.5.2a" \
"scikit-learn=0.18.1" \
"matplotlib=1.5.1"

echo -e "\n### Installing graphviz ###\n"

# installing or searching for graphviz in the anaconda channel failed a few times.
# Seemed to be an ittermittent anaconda channel server issue.

#conda search -c anaconda "graphviz"

conda install -y \
-c anaconda \
"graphviz=2.38.0=5"

# this link might work, just extract into miniconda folder and hope for the best
# - nope, missing libltdl.so.7
#cd ${BASE_DIR}/internals/thirdparty/miniconda
#p='https://anaconda.org/anaconda/graphviz/2.38.0/download/linux-64/graphviz-2.38.0-5.tar.bz2'
#curl -L ${p} -o graphviz-2.38.0-5.tar.bz2
#tar xjf graphviz-2.38.0-5.tar.bz2
# note: this will overwrite some libs, unsure if it will break environment

cd ${BASE_DIR}/internals

#source "${BASE_DIR}/internals/install/thirdparty_paths.sh"
#which dot
#dot -V
#exit


# anaconda graphviz installed after matplotlib 2.1.0 causes matplotlib to be downgraded
# (conda-forge superceded by anaconda channel)
# matplotlib: 2.1.1-py35_0          conda-forge --> 2.0.2-np113py35_0
# install matplotlib after graphviz to workaround hopefully? 
# - seems to somehow bring in an openssl lib requiring a newer GLIBC
# - disabling auto_update_conda seems to prevent unwanted openssl updates

# even if pinned, matplotlib and libgcc superceded by anaconda channel
# screw it, just going back to matplotlib 1.5.1 for now and backporting some rcParams into plotting scripts

#------------------------------------------------------------------------------------------------------
# cython executables in mpl 2.1.0 build are not binary compatible with environment, but not obvious until
# attempting to import matplotlib.pyplot
#mv ${CONDA_PATH}/lib/libcgraph.so.6.0.0 ${CONDA_PATH}/lib/libcgraph.so.6.0.0.bk
#mv ${CONDA_PATH}/lib/libgvc.so.6.0.0 ${CONDA_PATH}/lib/libgvc.so.6.0.0.bk

#echo -e "\n### Installing matplotlib ###\n"

#conda install -y \
#-c conda-forge \
#"matplotlib=2.1.0"

#mv ${CONDA_PATH}/lib/libcgraph.so.6.0.0.bk ${CONDA_PATH}/lib/libcgraph.so.6.0.0
#mv ${CONDA_PATH}/lib/libgvc.so.6.0.0.bk ${CONDA_PATH}/lib/libgvc.so.6.0.0
#------------------------------------------------------------------------------------------------------


echo -e "\n### Reinstalling libgcc ###\n"

# hack reinstall libgcc to workaround collision with libstdcxx-ng (brought in by matplotlib=2.1.0 and/or graphviz)
# --no-deps?
conda install -y -c defaults --force "libgcc=5.2.0"


# alternatively, just build it ourselves (requires cairo and pango)
#cd ${BASE_DIR}/thirdparty
#wget http://www.graphviz.org/pub/graphviz/stable/SOURCES/graphviz-2.40.1.tar.gz
#tar -xvf graphviz-2.40.1.tar.gz
#cd graphviz-2.40.1
#./configure
#make
#cd cmd/dot
#make dot_static


# - conda stuff interacting in weird ways with compile process, so remove from PATH
# from above: export PATH=${CONDA_PATH}/bin:$PATH
# so now just remove the first element in PATH
export PATH=${PATH#*:}

# bowtie2 more recent release
cd ${BASE_DIR}/internals/thirdparty
p='https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip'
if [ ! -f 'bowtie2-2.3.4.3-linux-x86_64.zip' ]; then
    echo "Downloading bowtie2"
    curl -L ${p} -o bowtie2-2.3.4.3-linux-x86_64.zip
fi
if [ ! -d bowtie2-2.3.4.3-linux-x86_64 ]; then
    unzip bowtie2-2.3.4.3-linux-x86_64.zip
fi
if [ ! -d 'bowtie2' ]; then
    mv bowtie2-2.3.4.3-linux-x86_64 bowtie2
fi

# ghostscript
cd ${BASE_DIR}/internals/thirdparty
p='https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs925/ghostscript-9.25-linux-x86_64.tgz'
if [ ! -f 'ghostscript-9.25.tar.gz' ]; then
    echo "Downloading ghostscript"
    curl -L ${p} -o ghostscript-9.25.tar.gz
fi
if [ ! -d ghostscript-9.25 ]; then
    tar -xf ghostscript-9.25.tar.gz
fi
if [ ! -d 'ghostscript' ]; then
    mv ghostscript-9.25-linux-x86_64 ghostscript
    cd ghostscript
    ln -s gs-925-linux-x86_64 gs 
fi


# pipeviewer                                                                                            
cd ${BASE_DIR}/internals/thirdparty
p='http://www.ivarch.com/programs/sources/pv-1.6.0.tar.gz'
if [ ! -f 'pv-1.6.0.tar.gz' ]; then
    echo "Downloading Pipe Viewer utility"
    #curl -L ${p} -o pv-1.6.0.tar.gz
    wget ${p} -O pv-1.6.0.tar.gz
fi
if [ ! -d pv-1.6.0 ]; then
    tar -xf pv-1.6.0.tar.gz
fi
if [ ! -d 'pipeviewer' ]; then
    cd pv-1.6.0
    ./configure --prefix=${BASE_DIR}/internals/thirdparty/pipeviewer
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

source "${BASE_DIR}/internals/paths/thirdparty_paths.sh"

which bowtie2-build
bowtie2-build --version
which bowtie2
bowtie2 --version
which bbmerge.sh

# bbmerge.sh --version exits with returncode 1, so just check string instead
set +e
result=$(bbmerge.sh --version 2>&1)
if [[ $result != *"BBMap version"* ]]; then
    echo "Error: could not detect BBmerge version. "
    echo "$result"
    exit 1
fi
set -e

which pv
pv --version
which python3.5
python3.5 --version
test_matplotlib
which dot
dot -V
which gs
gs --version



