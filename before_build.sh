#!/bin/bash
set -e -u -x

if [ "$(uname)" == "Linux" ] && [ $# -lt 1 ]; then
    cd /project
fi

git clone https://github.com/libarchive/libarchive.git
cd libarchive
cmake -DCMAKE_BUILD_TYPE=Release .
make -j4
if [ "$(uname)" == "Darwin" ]; then
    sudo make install
    export LDFLAGS="-L/usr/local/lib" 
    export CPPFLAGS="-I/usr/local/include -std=c++11"
else
    make install
fi

cd ..

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. -DLA_PREFIX="/usr/local"  
make -j4
cd ..
