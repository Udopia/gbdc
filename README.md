# GBDHash

* Efficiently calculate the GBD-Hash for a given DIMACS CNF file
* Uses libarchive to read compressed files, so it supports a wide range of compressed formats

## Standalone Application `gdbhash`

Run `make` to build. Requires presence of `libarchive` (e.g. for Debian/Ubuntu/Mint run `apt install libarchive-dev` first). 

## Python Module `gbdhashc`

This is mainly here to be integrated into GBD-Tools, because byte-wise processing of large files is ridiculously slow in Python. Run `make; make install` in order to setup the pyton module `gbdhashc`. Installations of GBD-Tools will then automatically start using this module to calculate gbd-hashes. This is orders of magnitude faster than the fallback implementation in Python. 