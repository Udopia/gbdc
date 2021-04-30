all: pymod

pymod: obj/%.o
	python3 setup.py build

install:
	python3 setup.py install
	
uninstall:
	python3 setup.py install --record files.tmp
	tr '\n' '\0' < files.tmp | xargs -0 sudo rm -f --
	rm files.tmp


