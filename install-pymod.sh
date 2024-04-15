#!/bin/bash

if [ $# -gt 0 ]; then
	j=$1
else
	j=1
fi

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$j
cd ..
pip install . --force-reinstall
rm -rf gbdc.egg-info
