# CNF to KIS transformer

- Source Context: cnf
- Target Context: kis

The CNF to KIS transformer translates SAT instances given in conjunctive normal form into instances of the k-Independent Set Problem.
The KIS instances are graphs given in DIMACS format with the following header "p kis $nodes $edges $k".
Here, $nodes denotes the number of nodes, $edges the number of edges in the graph, and $k the number of minimum independent sets that must exist in the graph for the instance to be satisfiable.
Our translation differs slightly from that presented in [1] in that instead of first translating the graph into a 3-SAT problem to introduce a triangle per clause, we directly introduce a clique per clause.

[1] Christos H. Papadimitriou. Computational complexity. 1994.