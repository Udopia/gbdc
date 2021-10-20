#mkdir -p build
#cd build
#cmake ..
#make
#cd ..
python3 setup.py build
python3 setup.py install --record uninstall.info --force
