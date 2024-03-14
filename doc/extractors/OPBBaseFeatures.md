# Base Features for OPB Instances

**Valid contexts**: `opb`

This feature extractor computes several statistical features for WCNF files.

## Feature Summary

| Category           | Feature name        | Description                                                                          |
|--------------------|---------------------|--------------------------------------------------------------------------------------|
| Basic statistics   | `constraints`       | Number of constraints                                                                |
|                    | `variables`         | Number of variables                                                                  |
| Constraint types   | `clauses`           | Number of clausal constraints                                                        |
|                    | `cards_ge`          | Number of cardinality upper bound (`>=`) constraints                                 |
|                    | `cards_eq`          | Number of cardinality equality (`=`) constraints                                     |
|                    | `pbs_ge`            | Number of pseudo-Boolean upper bound (`>=`) constraints                              |
|                    | `pbs_eq`            | Number of pseudo-Boolean equality (`=`) constraints                                  |
|                    | `assignments`       | Number of constraints admitting only a single assignment for the variables appearing |
|                    | `trivially_unsat`   | Number of trivially unsatisfiable constraints                                        |
| Objective function | `obj_terms`         | Number of terms in the objective function                                            |
|                    | `obj_max_val`       | Maximum value that the objective can take                                            |
|                    | `obj_min_val`       | Minimum value that the objective can take                                            |
|                    | `obj_coeffs_<dist>` | Distribution of objective coefficients                                               |

Distributions features are `mean`, `variance`, `min`, `max` and `entropy`.
