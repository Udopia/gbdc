mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
pip install . --user --force-reinstall
rm -rf gbdc.egg-info