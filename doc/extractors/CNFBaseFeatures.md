# Base Features for CNF Instances

**Valid contexts**: `cnf`

This feature extractor computes several statistical features which can roughly be grouped into the four categories _clause type_, _horn proximity_, _polarity balance_, and _graph based_.

## Feature Summary

| Category              | Feature Name           | Description                                                                     |
|-----------------------|------------------------|---------------------------------------------------------------------------------|
| Basic statistics      | `clauses`              | Number of clauses                                                               |
|                       | `variables`            | Number of variables                                                             |
|                       | `bytes`                | Number of bytes/characters in text file                                          |
|                       | `ccs`                  | Number of connected components                                                  |
| Clause types          | `cls<n>` for `<n>` 1-9 | Number of clauses with `<n>` literals                                           |
|                       | `cls10p`               | Number of clauses with 10 or more literals                                      |
|                       | `horn`                 | Number of horn clauses                                                          |
|                       | `invhorn`              | Number of inverted horn clauses                                                 |
|                       | `positive`             | Number of positive clauses                                                      |
|                       | `negative`             | Number of negative clauses                                                      |
| Horn proximity        | `hornvars_<dist>`      | Distribution of variable occurrence in horn clauses                             |
|                       | `invhornvars_<dist>`   | Distribution of variable occurrence in inverted horn clauses                    |
| Polarity balance      | `balancecls_<dist>`    | Distributions of fractions of number of positive/negative literals in clauses   |
|                       | `balancevars_<dist>`   | Distributions of fractions of number of positive/negative literals of variables |
| Variable-clause graph | `vcg_vdegree_<dist>`   | Variable degree distribution                                                    |
|                       | `vcg_cdegree_<dist`    | Clause degree distribution                                                      |
| Variable graph        | `vg_degree_<dist>`     | Degree distribution                                                             |
| Clause graph          | `cg_degree_<dist>`     | Degree distribution                                                             |

Distributions features are `mean`, `variance`, `min`, `max` and `entropy`.
