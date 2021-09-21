# CNF Tools

Standalone Application and Python Accelerator Module

* Efficiently calculate the GBD-Hash for a given DIMACS CNF file
* Can read a variety of compressed formats (using libarchive)
* Normalized output of DIMACS CNF 

## Standalone Application `cnftools`

Run `make` to build. Requires presence of `libarchive` (e.g. for Debian/Ubuntu/Mint run `apt install libarchive-dev` first). 

## Python Accelerator Module `gbdc`

Setup Accelerator Module for GBD Tools

Run `make; make install` in order to install the python accelerator module `gbdc`. 
GBD Tools uses the `gdbc` module if it is installed. 
