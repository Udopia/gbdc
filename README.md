# GBDC: Global Benchmark Database, C++ Extension Module

[![Linux Build-Test](https://github.com/Udopia/gbdc/actions/workflows/linux_build_test.yml/badge.svg?branch=master)](https://github.com/Udopia/gbdc/actions/workflows/linux_build_test.yml)
[![MacOS Build-Test](https://github.com/Udopia/gbdc/actions/workflows/macos_build_test.yml/badge.svg?branch=master)](https://github.com/Udopia/gbdc/actions/workflows/macos_build_test.yml)
[![Windows Build-Test](https://github.com/Udopia/gbdc/actions/workflows/windows_build_test.yml/badge.svg?branch=master)](https://github.com/Udopia/gbdc/actions/workflows/windows_build_test.yml)

[GBDC](https://github.com/Udopia/gbdc) provides efficient implementations of functions for benchmark instance identification, instance feature extraction and instance transformation.
GBDC provides a command-line tool as well as the Python package `gbdc`.
The Python package `gbdc` is used by [Global Benchmark Database](https://github.com/Udopia/gbd).

## Supported Domains and Features

GBDC provides benchmark instance identifiers, feature extractors, and instance transformers for several problem domains, including propositional satisfiability (SAT) and optimization (MaxSAT), as well as Pseudo-Boolean Optimization (PBO).
Inluded are the following implementations.
A documentation on the supported domains, feature extractors, and instance transformers can be found in the [documentation](https://udopia.github.io/gbdc/doc/Index.html).

## Build Dependencies and Installation Instructions

* GBDC uses `libarchive` for reading from a large variety of compressed formats (in some systems provided by the package `libarchive-dev`).
* Some GBDC functions use an [IPASIR](https://github.com/biotomas/ipasir) SAT Solver. GBDC's build-system pulls the external SAT Solver [CaDiCaL](http://fmv.jku.at/cadical/) by A. Biere (MIT licensed).

<!-- #### Shipped Dependencies

* A copy of the command-line argument parser by P. S. Kumar [`argparse.h`](https://github.com/p-ranav/argparse) (MIT licensed) resides in the `lib` folder.

* A copy of the [MD5 hash](https://github.com/CommanderBubble/MD5) implementation by M. Lloyd (MIT licensed) resides in the `lib` folder. -->

```bash
# Build command-line tool and dependencies for the Python module
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
# Install Python module
python3 setup.py build
python3 setup.py install --user --record uninstall.info
```

<!-- ## Publications

* Gate feature extraction uses our gate recognition algorithm which is described in the following publications:

    * [*Recognition of Nested Gates in CNF Formulas* (SAT 2015, Iser et al.)](https://rdcu.be/czCr1)

    * [*Recognition and Exploitation of Gate Structure in SAT Solving* (2020, Iser)](https://d-nb.info/1209199122/34)

* The Python module `gbdc` is used in our project [GBD Benchmark Database](https://github.com/Udopia/gbd)

    * [*Collaborative Management of Benchmark Instances and their Attributes* (2020, Iser et al.)](https://arxiv.org/pdf/2009.02995.pdf) -->
