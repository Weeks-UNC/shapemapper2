#!/bin/bash

#set -e #exit on error
#set -x #Run script in debug mode


eval "$(conda shell.bash hook)" #Shell configuration work around so conda activate works in slurm 
# Finds absolute path of shapemapper directory

# Solution found on stack overflow to find directory name if its a slurm job
if [ -n "$SLURM_JOB_ID" ];  then
   # check the original location through scontrol and $SLURM_JOB_ID
   BASE_DIR=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
   BASE_DIR=$( dirname $( dirname $( dirname $BASE_DIR )))
else
   # otherwise: started with bash. Get the real location.
   BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"
   #SCRIPT_PATH=$(realpath $0)
fi

## Writes output to a log file
##exec &> >(tee "$BASE_DIR/internals/build_thirdparty.log")

unset PYTHONPATH 
CONDA_PATH=${BASE_DIR}/internals/thirdparty/miniconda
export PATH=${CONDA_PATH}/bin:$PATH

echo "++++++++++++++++++++++"
echo "Conda path: $CONDA_PATH"
wconda=$(which conda)
echo "Which conda: $wconda"
echo "BASE DIR: ${BASE_DIR}"
echo "++++++++++++++++++++++"

echo "Creating internal conda environment"
conda activate $CONDA_PATH
conda env create --prefix=${BASE_DIR}/internals/thirdparty/miniconda/envs/shapemapper_make -f $BASE_DIR/internals/thirdparty_helper/shapemapper_compile_env.yml


conda deactivate
echo "Activating internal conda environment"
conda activate ${BASE_DIR}/internals/thirdparty/miniconda/envs/shapemapper_make

cd ${BASE_DIR}/internals/thirdparty


#Installs bowtie2
cd ${BASE_DIR}/internals/thirdparty
p='https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip'
if [ ! -f 'bowtie2-2.3.4.3-linux-x86_64' ];
then
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

#Installs pipeviewer
cd ${BASE_DIR}/internals/thirdparty
p='http://www.ivarch.com/programs/sources/pv-1.6.20.tar.gz'
if [ ! -f 'pv-1.6.20.tar.gz' ]; then
   echo "Downloading pipe viewer utility"
   #wget ${p} -o pv-1.6.20.tar.gz
   curl -L ${p} -o pv-1.6.20.tar.gz
fi
if [ ! -d pv-1.6.20 ]; then
   tar xvf pv-1.6.20.tar.gz
fi
if [ ! -d 'pipeviewer' ]; then
   cd pv-1.6.20
   ./configure --prefix=${BASE_DIR}/internals/thirdparty/pipeviewer
   make
   make install
fi


#build bbmerge
export JAVA_HOME=${CONDA_PATH}/envs/shapemapper_make/conda-meta/java-jdk-8.0.112-1.json
cd ${BASE_DIR}/internals/thirdparty/miniconda/envs/shapemapper_make/opt/bbmap-37.78/jni/
make -f makefile.linux


source ${BASE_DIR}/internals/paths/thirdparty_paths.sh
echo "------------Running checks-------------"
which bowtie2-build
bowtie2-build --version
which bowtie2
bowtie2 --version
which bbmerge.sh
result=$(bbmerge.sh --version 2>&1)
set +e
if [[ $result != *"BBMap version"* ]]; then
    echo "Error: could not detect BBmerge version. "
    echo "$result"
    exit 1
fi
set -e

which pv
pv --version
which python3.9
python3.9 --version
which gs
gs --version
