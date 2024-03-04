[![Linux Build-Test](https://github.com/Udopia/gbdc/actions/workflows/linux_build_test.yml/badge.svg?branch=master)](https://github.com/Udopia/gbdc/actions/workflows/linux_build_test.yml)
[![MacOS Build-Test](https://github.com/Udopia/gbdc/actions/workflows/macos_build_test.yml/badge.svg?branch=master)](https://github.com/Udopia/gbdc/actions/workflows/macos_build_test.yml)
[![Windows Build-Test](https://github.com/Udopia/gbdc/actions/workflows/windows_build_test.yml/badge.svg?branch=master)](https://github.com/Udopia/gbdc/actions/workflows/windows_build_test.yml)

# Global Benchmark Database, C++ Extension Module

This project provides the Python module `gbdc`, which offers efficient implementations of functions for benchmark instance identification, feature extraction and problem transformation in C++.
These are used by the initialization functions of the database [GBD](https://github.com/Udopia/gbd).
The functions are also available via the command-line tool `gbdc`, which can also be created from this repository.

## Supported Domains and Features

GBDC provides benchmark instance identifiers, feature extractors, and instance transformers for several problem domains, including those of propositional satisfiability (SAT) and optimization (MaxSAT), as well as Pseudo-Boolean Optimization (PBO).
Inluded are the following implementations.

* Implementations of `GBDHash` for domain-specific instance identifier calculation.
* Implementations of `ISOHash` for domain-specific isomorphism class identification.
* Implementations of several domain-specific feature extractors.
* Implementations of instance normalizers and problem transformers.

More information on the supported domains can be found in the [documentation of domains](doc/Domains.md).
More information on the supported feature extractors can be found in the [documentation of feature extractors](doc/Extractors.md).
More information on the supported problem transformers can be found in the [documentation of problem transformers](doc/Extractors.md).

## Build Dependencies and Installation Instructions

* GBDC uses `libarchive` for reading from a large variety of compressed formats (in some systems provided by the package `libarchive-dev`).
* Some GBDC functions use an [IPASIR](https://github.com/biotomas/ipasir) SAT Solver. GBDC's build-system pulls the external SAT Solver [CaDiCaL](http://fmv.jku.at/cadical/) by A. Biere (MIT licensed).

<!-- #### Shipped Dependencies

* A copy of the command-line argument parser by P. S. Kumar [`argparse.h`](https://github.com/p-ranav/argparse) (MIT licensed) resides in the `lib` folder.

* A copy of the [MD5 hash](https://github.com/CommanderBubble/MD5) implementation by M. Lloyd (MIT licensed) resides in the `lib` folder. -->

```bash
# Build command-line tool
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
