# Isohash

Isohash over-approximates the equivalence class of isomorphic instances in their context.
An implementation for Isohash exists for several contexts.

## Isohash for CNF Instances

- context: cnf

This feature extractor computes exactly one feature for each instance which is the hashsum of the sorted degree sequence of the normalize literal interaction graph. TODO: elaborate.
The feature over-approximates the equivalence class of isomorphic instances and can be used to reduce redundancy in experimentation.

## Isohash for WCNF Instances

- context: wcnf

