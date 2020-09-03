# GBDHash

* Efficiently calculate the GBD-Hash for a given DIMACS CNF file
* Uses libarchive to read compressed files, so it supports a wide range of compressed formats

This is mainly here to be integrated into GBD, because byte-wise processing of large files is ridiculously slow in Python. If you want to use it as a standalone application then check it out, run make, and give it a go.