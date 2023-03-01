#!/usr/bin/bash
set -e # exit on first error (if any)
eval "$(conda shell.bash hook)"
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "-------"
echo "Base dir: $BASE_DIR"
echo "-------"
CONDA_PATH=$BASE_DIR/../../internals/thirdparty/miniconda
export PATH=${CONDA_PATH}/bin:$PATH


echo "-------checking conda location-------" 
which conda

echo "-------activating conda location-------" 
conda activate $BASE_DIR/../../internals/thirdparty/miniconda/envs/shapemapper_make
conda env list
which python
which make
which cmake
python --version

echo "------building binaries-------" 
cd ../..
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
make
