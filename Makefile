all: gbdhash pymod

gbdhash: Main.cc StreamBuffer.h md5/md5.cpp md5/md5.h md5/md5_loc.h
	g++ -o gbdhash Main.cc md5/md5.cpp -larchive

pymod: Main.cc StreamBuffer.h md5/md5.cpp md5/md5.h md5/md5_loc.h
	python3 setup.py build

install:
	python3 setup.py install
	
uninstall:
	python3 setup.py install --record files.tmp
	tr '\n' '\0' < files.tmp | xargs -0 sudo rm -f --
	rm files.tmp

clean:
	rm gbdhash
	rm -rf build

