# Isohash

`isohash` over-approximates the equivalence class of isomorphic instances in their context.
An implementation for Isohash exists for several contexts.

## Isohash for CNF Instances

**Valid contexts**: `cnf`

This feature extractor computes exactly one feature (`isohash`) for each instance which is the hashsum of the sorted degree sequence of the normalized literal interaction graph.
The feature over-approximates the equivalence class of isomorphic instances and can be used to reduce redundancy in experimentation.

In detail, the feature is computed by generating a literal interaction graph based on the clauses in the instance.
For variables where the negative literal degree is higher than the positive one, the degrees are flipped.
Next, a variable order is determined based on the positive literal degree, using the negative degree as a tie breaker.
Lastly, the degree sequence is hashed with an MD5 hash, following the determined variable order.
The final hashed string looks something like `<v1 pos deg> <v1 neg deg> <v2 pos deg> <v2 neg deg> ...` with `v1`, `v2` being the variables in the determined order and `pos/neg deg` being their positive/negative literal degree in the graph.

## Isohash for WCNF Instances

**Valid contexts**: `wcnf`

This feature extractor computes exactly one feature (`isohash`) for each instance.
The feature over-approximates the equivalence class of isomorphic instances and can be used to reduce redundancy in experimentation.
The hash is built based on the literal interaction graphs, one taking only the hard clauses into account and one taking all clauses into account.
For the graph taking all clauses into account, hard clauses are treated with edge weight 1, while soft clauses use their weight as the edge weight.
Finally, the two degree sequences are sorted as for CNF files and hashed while being joined together by the literal `"softs "`.
