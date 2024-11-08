#!/bin/bash

install_pybind_devel_from_source() {
    git clone https://github.com/pybind/pybind11.git
    cd pybind11
    mkdir build
    cd build
    cmake ..
    make -j6 install
}

install_libarchive_devel_from_source() {
    http://http.debian.net/debian/pool/main/liba/libarchive/libarchive_3.2.2.orig.tar.gz
    tar -zxvf libarchive_3.2.2.orig.tar.gz
    cd libarchive_3.2.2.orig
    ./configure
    make
    sudo make install
}

if command -v dnf > /dev/null; then
    dnf install -y libarchive-devel
    dnf install -y cmake
    dnf groupinstall -y "Development Tools"
    install_pybind_devel_from_source
elif command -v yum > /dev/null; then
    #yum install -y libarchive-devel
    yum install -y cmake
    yum groupinstall -y "Development Tools"
    install_pybind_devel_from_source
    install_libarchive_devel_from_source
else
    apt update
    apt install -y libarchive-dev
    apt install -y cmake
    apt install -y build-essential
    apt install -y pybind11-dev
fi