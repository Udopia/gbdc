#!/bin/bash

install_pybind_devel_from_source() {
    git clone https://github.com/pybind/pybind11.git
    cd pybind11
    mkdir build
    cd build
    cmake ..
    make -j6 install
}

if command -v dnf > /dev/null; then
    dnf install -y libarchive-devel
    dnf install -y cmake
    dnf groupinstall -y "Development Tools"
    install_pybind_devel_from_source
elif command -v yum > /dev/null; then
    yum install -y libarchive-devel
    yum install -y cmake
    yum groupinstall -y "Development Tools"
    install_pybind_devel_from_source
else
    apt update
    apt install -y libarchive-dev
    apt install -y cmake
    apt install -y build-essential
    apt install -y pybind11-dev
fi