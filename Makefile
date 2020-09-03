all: test

test: Main.cc StreamBuffer2.h md5/md5.cpp md5/md5.h md5/md5_loc.h
	g++ -o gbdhash Main.cc md5/md5.cpp -larchive

clean:
	rm gbdhash

