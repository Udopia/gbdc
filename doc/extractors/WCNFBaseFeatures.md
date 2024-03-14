# Base Features for WCNF Instances

**Valid contexts**: `wcnf`

This feature extractor computes several statistical features for WCNF files.
The features related to the hard clauses are the same ones computed in the [CNF base feature extractor](CNFBaseFeatures.md).
Additional features related to the soft clauses are added.
Features related to the hard clauses are prefixed with `h_` and features related to the soft clauses with `s_`.

## Feature Summary

| Category                   | Feature name           | Description                                                                                                   |
|----------------------------|------------------------|---------------------------------------------------------------------------------------------------------------|
| Basic statistics           | `h_clauses`            | Number of hard clauses                                                                                        |
|                            | `variables`            | Number of variables                                                                                           |
|                            | `s_clauses`            | Number of soft clauses                                                                                        |
|                            | `s_weight_sum`         | Sum of weights                                                                                                |
| Hard clause types          | `h_cls<n>`             | Number of hard clauses of length `<n>` (1-9)                                                                  |
|                            | `h_cls10p`             | Number of hard clauses with at least 10 literals                                                              |
|                            | `h_horn`               | Number of hard horn clauses                                                                                   |
|                            | `h_invhorn`            | Number of hard inverted horn clauses                                                                          |
|                            | `h_positive`           | Number of hard positive clauses                                                                               |
|                            | `h_negative`           | Number of hard negative clauses                                                                               |
| Hard horn proximity        | `h_hornvars_<dist>`    | Distribution of variable occurrences in hard horn clauses                                                     |
|                            | `h_invhornvars_<dist>` | Distribution of variable occurrences in hard inverse horn clauses                                             |
| Hard polarity balance      | `h_balancecls_<dist>`  | Distributions of fractions of number of positive/negative literals in hard clauses                            |
|                            | `h_balancevars_<dist>` | Distributions of fractions of number of positive/negative literals of variables only considering hard clauses |
| Hard variable-clause graph | `h_vcg_vdegree_<dist>` | Variable degree distribution                                                                                  |
|                            | `h_vcg_cdegree_<dist>` | Clause degree distribution                                                                                    |
| Hard variable graph        | `h_vg_degree_<dist>`   | Degree distribution                                                                                           |
| Hard clause graph          | `h_cg_degree_<dist>`   | Degree distribution                                                                                           |
| Soft clause types          | `s_cls<n>`             | Number of soft clauses of length `<n>` (1-9)                                                                  |
|                            | `s_cls10p`             | Number of soft clauses with at least 10 literals                                                              |
| Weight distribution        | `s_weight_<dist>`      | Distribution of weights                                                                                       |

Distributions features are `mean`, `variance`, `min`, `max` and `entropy`.
