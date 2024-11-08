#!/bin/bash

install_pybind_devel_from_source() {
    git clone https://github.com/pybind/pybind11.git
    cd pybind11
    mkdir build
    cd build
    cmake ..
    make -j6 install
    cd ../..
}

install_libarchive_devel_from_source() {
    #git clone https://github.com/libarchive/libarchive.git
    #cd libarchive
    #/bin/sh build/autogen.sh
    curl -L https://libarchive.org/downloads/libarchive-3.7.7.tar.gz -o libarchive-3.7.7.tar.gz
    tar -xvf libarchive-3.7.7.tar.gz
    cd libarchive-3.7.7
    ./configure
    make
    make check
    make install
    cd ..
}

if command -v dnf > /dev/null; then
    dnf update
    dnf install -y cmake
    dnf groupinstall -y "Development Tools"
    install_pybind_devel_from_source
    install_libarchive_devel_from_source
elif command -v yum > /dev/null; then
    yum update
    yum install -y cmake
    yum groupinstall -y "Development Tools"
    install_pybind_devel_from_source
    install_libarchive_devel_from_source
elif command -v apt > /dev/null; then
    apt update
    apt install -y libarchive-dev
    apt install -y cmake
    apt install -y build-essential
    apt install -y pybind11-dev
elif command -v apk > /dev/null; then
    apk update
    apk install -y libarchive-dev
    apk install -y cmake
    apk install -y build-base
    install_pybind_devel_from_source
else
    echo "No package manager found"
    exit 1
fi