# Gate Features for CNF Instances

**Valid contexts**: `cnf`

Executes the gate extraction algorithm described in [[1]](https://nbn-resolving.org/urn:nbn:de:101:1-2020042904595660732648) and extracts features that indicate the number, type, and position of gates in the extracted hierarchical gate structure.

## Feature Summary

"levels_mean", "levels_variance", "levels_min", "levels_max", "levels_entropy"
"levels_none_mean", "levels_none_variance", "levels_none_min", "levels_none_max", "levels_none_entropy"
"levels_generic_mean", "levels_generic_variance", "levels_generic_min", "levels_generic_max", "levels_generic_entropy"
"levels_mono_mean", "levels_mono_variance", "levels_mono_min", "levels_mono_max", "levels_mono_entropy"
"levels_and_mean", "levels_and_variance", "levels_and_min", "levels_and_max", "levels_and_entropy"
"levels_or_mean", "levels_or_variance", "levels_or_min", "levels_or_max", "levels_or_entropy"
"levels_triv_mean", "levels_triv_variance", "levels_triv_min", "levels_triv_max", "levels_triv_entropy"
"levels_equiv_mean", "levels_equiv_variance", "levels_equiv_min", "levels_equiv_max", "levels_equiv_entropy"
"levels_full_mean", "levels_full_variance", "levels_full_min", "levels_full_max", "levels_full_entropy"

| Category              | Feature Name           | Description                                                                     |
|-----------------------|------------------------|---------------------------------------------------------------------------------|
| Basic statistics      | `n_vars`               | Number of variables                                                             |
|                       | `n_gates`              | Number of gates                                                                 |
|                       | `n_roots`              | Number of roots / output gates                                                  |
|                       | `n_none`               | Number of input variables                                                       |
| Gate types            | `n_mono`               | Number of monotonically nested gates                                            |
|                       | `n_and`                | Number of AND gates                                                             |
|                       | `n_or`                 | Number of OR gates                                                              |
|                       | `n_triv`               | Number of unary equivalences                                                    |
|                       | `n_equiv`              | Number of binary equivalence and XOR gates                                      |
|                       | `n_full`               | Number of gates detected by the 2^n all-different clauses criterion             |
|                       | `n_generic`            | Number of gates detected with the semantic criterion by a SAT solver            |
| Position of Gates     | `levels_<dist>`        | Levels of gates in the reconstructed DAG                                        |
|                       | `levels_none_<dist>`   | Levels of input variables                                                       |
|                       | `levels_mono_<dist>`   | Levels of monotonically nested gates                                            |
|                       | `levels_and_<dist>`    | Levels of AND gates                                                             |
|                       | `levels_or_<dist>`     | Levels of OR gates                                                              |
|                       | `levels_triv_<dist>`   | Levels of unary equivalences                                                    |
|                       | `levels_equiv_<dist>`  | Levels of binary equivalences                                                   |
|                       | `levels_full_<dist>`   | Levels of gates detected by the 2^n all-different clauses criterion             |
|                       | `levels_generic_<dist>`| Levels of gates detected with the semantic criterion by a SAT solver            |

Distributions features are `mean`, `variance`, `min`, `max` and `entropy`.

## References

- [1] Markus Iser: _Recognition and Exploitation of Gate Structure in SAT Solving_, 2020.
