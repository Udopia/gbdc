#!/bin/bash

dnf install libarchive-devel
dnf install cmake
dnf groupinstall "Development Tools"
git clone https://github.com/pybind/pybind11.git
cd pybind11
mkdir build
cd build
cmake ..
make install