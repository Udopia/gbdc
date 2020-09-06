all: gbdhash pymod

gbdhash: Main.cc StreamBuffer2.h md5/md5.cpp md5/md5.h md5/md5_loc.h
	g++ -o gbdhash Main.cc md5/md5.cpp -larchive

pymod: Main.cc StreamBuffer2.h md5/md5.cpp md5/md5.h md5/md5_loc.h
	python3 setup.py build

install:
	python3 setup.py install

clean:
	rm gbdhash
	rm -rf build

