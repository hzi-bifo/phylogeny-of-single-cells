#!/usr/bin/env bash
set -euo pipefail

git clone --recurse-submodules -j 2 https://github.com/amkozlov/raxml-ng
cd raxml-ng
git checkout 816f9d1395871443607a727d1908ad96de993d2e
mkdir build && cd build
cmake ..
make
install -d ${CONDA_PREFIX}/bin
install ../bin/raxml-ng ${CONDA_PREFIX}/bin
cd ../..
rm -rf raxml-ng