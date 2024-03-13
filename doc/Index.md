# Supported Benchmark Problem Domains and Formats

GBDC reads benchmark instances from text file based formats.
Those text files are assigned to the respective domain by the associated extension.
Text files can also be packed, in that case the extension can be augmented by .xz, .lzma, .bz2, or .gz.

## Propositional Satisfiability (SAT)

* Context: cnf
* Format: DIMACS CNF (todo: add link or reference)
* Extensions: .cnf, .cnf.xz, .cnf.bz2, .cnf.lzma, .cnf.gz
* Identifier:
  * no line-breaks or header
  * space separated literals
  * zero separated clauses
  * md5 hash
* Feature Extractors:
  * [Isohash](extractors/Isohash.md)
  * [Base Features](extractors/CNFBaseFeatures.md)
  * [Gate Features](extractors/CNFGateFeatures.md)

## Propositional Optimization (MaxSAT)

* Context: wcnf
* Format: DIMACS WCNF (todo: add link or reference)
* Extensions: .wcnf, .wcnf.xz, .wcnf.bz2, .wcnf.lzma, .wcnf.gz
* Identifier: [todo]
* Feature Extractors:
  * [Isohash](extractors/Isohash.md)
  * [Base Features](extractors/WCNFBaseFeatures.md)

## Pseudo-Boolean Optimization (PBO)

* Context: opb
* Format: OPB Format (todo: add link or reference)
* Extension: .opb, .opb.xz, .opb.bz2, .opb.lzma, .opb.gz
* Identifier: [todo]
  * [Base Features](extractors/OPBBaseFeatures.md)
