#!/bin/bash

if command -v dnf > /dev/null; then
    dnf install -y libarchive-devel
    dnf install -y cmake
    dnf groupinstall -y "Development Tools"
    git clone https://github.com/pybind/pybind11.git
    cd pybind11
    mkdir build
    cd build
    cmake ..
    make -j6 install
else
    sudo apt update
    sudo apt install -y libarchive-dev
    sudo apt install -y cmake
    sudo apt install -y build-essential
    sudo apt install -y pybind11-dev
fi