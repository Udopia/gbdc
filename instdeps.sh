#!/bin/bash

dnf install -y libarchive-devel
dnf install -y cmake
dnf groupinstall -y "Development Tools"
# cd /project
git clone https://github.com/pybind/pybind11.git
cd pybind11
mkdir build
cd build
cmake ..
make -j6 install