# Supported Benchmark Problem Domains and Formats

GBDC reads benchmark instances from text file based formats.
Those text files are assigned to the respective domain by the associated extension.
Text files can also be packed, in that case the extension can be augmented by .xz, .lzma, .bz2, or .gz.

## Propositional Satisfiability (SAT)

* Format: DIMACS CNF
* Extension: .cnf
* Identifier:
  * no line-breaks or header
  * space separated literals
  * zero separated clauses
  * md5 hash

## Propositional Optimization (MaxSAT)

* Format: WCNF
* Extension: .wcnf

## Pseudo-Boolean Optimization (PBO)

* Format: OPB
* Extension: .opb
