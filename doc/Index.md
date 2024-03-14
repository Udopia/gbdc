# Supported Benchmark Problem Domains and Formats

GBDC reads benchmark instances from text file based formats.
Those text files are assigned to the respective domain by the associated extension.
Text files can also be packed, in that case the extension can be augmented by `.xz`, `.lzma`, `.bz2`, or `.gz`.

## Propositional Satisfiability (SAT)

* Context: `cnf`
* Format: [DIMACS CNF](http://www.satcompetition.org/2011/format-benchmarks2011.html)
* Extensions: `.cnf`, `.cnf.xz`, `.cnf.bz2`, `.cnf.lzma`, `.cnf.gz`
* Identifier:
  * no line breaks, header, or comments
  * single space separated literals
  * zero separated clauses
  * [MD5 hash](https://en.wikipedia.org/wiki/MD5)
* Feature Extractors:
  * [Isohash](extractors/Isohash.md)
  * [Base Features](extractors/CNFBaseFeatures.md)
  * [Gate Features](extractors/CNFGateFeatures.md)
* Transformers:
  * [CNF2KIS](transformers/CNF2KIS.md)
  * [Sanitizer](transformers/Sanitizer.md)

## Propositional Optimization (MaxSAT)

* Context: `wcnf`
* Format: [DIMACS WCNF](https://maxsat-evaluations.github.io/2024/rules.html#input) (supports both pre-2022 and new format)
* Extensions: `.wcnf`, `.wcnf.xz`, `.wcnf.bz2`, `.wcnf.lzma`, `.wcnf.gz`
* Identifier:
  * convert to new WCNF format
  * no line breaks, header, or comments
  * space separated literals
  * [MD5 hash](https://en.wikipedia.org/wiki/MD5)
* Feature Extractors:
  * [Isohash](extractors/Isohash.md)
  * [Base Features](extractors/WCNFBaseFeatures.md)

## Pseudo-Boolean Optimization (PBO)

* Context: `opb`
* Format: [OPB Format](https://www.cril.univ-artois.fr/PB07/solver_req.html)
* Extension: `.opb`, `.opb.xz`, `.opb.bz2`, `.opb.lzma`, `.opb.gz`
* Identifier:
  * no line breaks, header, or comments
  * remove unnecessary spaces
  * [MD5 hash](https://en.wikipedia.org/wiki/MD5)
* Feature Extractors:
  * [Base Features](extractors/OPBBaseFeatures.md)
