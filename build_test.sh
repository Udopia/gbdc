#!/bin/bash

if [ "$1" != "Debug" ] && [ "$1" != "Release" ]; then
    echo "Missing build type: Please specify 'Debug' or 'Release' as the first argument."
    exit 1
fi

build_type=$1

if [[ $# -gt 1 ]]; then
    j=$2
else
    j=1
fi

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=$build_type ..
make -j$j
ctest .
cd ..
