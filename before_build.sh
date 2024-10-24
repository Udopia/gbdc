#!/bin/bash
set -e -u -x

if [ "$(uname)" == "Linux" ]; then
    dnf install libarchive-devel
    dnf install cmake
    dnf groupinstall "Development Tools"
    git clone https://github.com/pybind/pybind11.git
    cd pybind11
    mkdir build
    cd build
    cmake ..
    make install
else
    brew install libarchive pybind11
fi