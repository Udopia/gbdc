# CNF Tools

CNF Tools provide the command-line program `cnftools` and the python module `gdbc`. All provided functionality is usually accessible from both the command-line but also from python via `gdbc`.

## Programming Language
- C++

## Tools

Input to all tools is a SAT instances given as a DIMACS CNF file which can be given in a variety of compressed formats (supported by libarchive). The following tools are provided:

* GBDHash:
> Calculates the identifier for the given instance which is used in [GBD Tools](https://pypi.org/project/gbd-tools/) for data organization. GBD Tools themselves use the provided python module `gdbc` if installed (with priority over its own fallback implementation in Python).
* Feature Extractors:
    * Base Features: The features cover degree distributions of well-known graph representations of a given instance and many more (see code for details).

    * Gate Features: The features cover gate distribuations over levels of the (potentially recoverable) hierarchical gate strucuture of an instance (see code for details).

* Problem Transformers:
    * ~~Sanitizer for DIMACS CNF (correct header, remove comments and extra whitespace, remove redundant literals in clause, delete tautological clauses)~~ (Tool not ready atm)

    * Transformation to Independent Set Problem


## Dependencies

* CNF Tools use `libarchive` to read instances from a large variety of compressed formats. Debian-based systems (incl. Ubuntu or Mint) provide this library in the package `libarchive-dev`.

* The program `cnftools` uses [`argparse.h` by Pranav Srinivas Kumar](https://github.com/p-ranav/argparse) to parse Command-Line Arguments. This is under MIT licence and resides in CNF Tools `lib` folder. It requires C++17 (and gcc 8+) to build.

* The tool `gbdhash` uses the [MD5 Hash implementation of Micheal Flyod](https://github.com/CommanderBubble/MD5). This is under MIT licence and resides in CNF Tools `lib` folder.

* Gate feature extraction uses a SAT Solver via the [IPSIR Interface](https://github.com/biotomas/ipasir). For this application, CNF Tools use the SAT Solver [CaDiCaL (by Armin Biere)](http://fmv.jku.at/cadical/) which underlies the MIT License.


## Publications

* Gate feature extraction uses our gate recognition algorithm which is described in the following publications:

    * [*Recognition of Nested Gates in CNF Formulas* (SAT 2015, Iser et al.)](https://rdcu.be/czCr1)

    * [*Recognition and Exploitation of Gate Structure in SAT Solving* (2020, Iser)](https://d-nb.info/1209199122/34)

* The Python module `gdbc` is used in our project [GBD Benchmark Database](https://github.com/Udopia/gbd)

    * [*Collaborative Management of Benchmark Instances and their Attributes* (2020, Iser et al.)](https://arxiv.org/pdf/2009.02995.pdf)


## Build and Install

### 1. Build `cnftools`:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

### 2. Install `gbdc`

    python3 setup.py build
    python3 setup.py install --record uninstall.info







