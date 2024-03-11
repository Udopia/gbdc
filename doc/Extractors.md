# Feature Extractors

## CNF Context

### Isohash for CNF Instances

This feature extractor computes exactly one feature for each instance which is the hashsum of the sorted degree sequence of the normalize literal interaction graph. TODO: elaborate. 
The feature over-approximates the equivalence class of isomorphic instances and can be used to reduce redundancy in experimentation.

### Base Features for CNF Instances

This feature extractor computes several statististical features which can roughly be grouped into the four categories \textsf{clause type}, \textsf{horn proximity}, \textsf{polarity balance}, and \textsf{graph based}.

- Clause type: includes clause counts of clauses for different fixed sizes, horn, inverted horn, as well as positive and negative clauses. The category encompasses the following features: 'clauses', 'variables', 'cls1', 'cls2', 'cls3', 'cls4', 'cls5', 'cls6', 'cls7',  'cls8', 'cls9', 'cls10p', 'horn', 'invhorn', 'positive', 'negative'

- Horn proximity: includes features describing the distribution of variable occurence counts in horn and inverted horn clauses. The category encompasses the following features: 'hornvars_mean', 'hornvars_variance', 'hornvars_min', 'hornvars_max', 'hornvars_entropy', 'invhornvars_mean', 'invhornvars_variance', 'invhornvars_min', 'invhornvars_max', 'invhornvars_entropy'

- Polarity balance: includes features describing the distribution of fractions of number of positive and negative literals in clauses and fractions of number of positive and negative literals of variables. The category encompasses the following features: 'balancecls_mean', 'balancecls_variance', 'balancecls_min', 'balancecls_max', 'balancecls_entropy', 'balancevars_mean', 'balancevars_variance', 'balancevars_min', 'balancevars_max', 'balancevars_entropy'

- Graph based: includes features describing the distribution of variable and clause degrees in the variable clause graph, the clause graph, and the variable graph. The category encompasses the following features:

  - based on the Variable-Clause Graph: 'vcg_vdegree_mean', 'vcg_vdegree_variance', 'vcg_vdegree_min', 'vcg_vdegree_max', 'vcg_vdegree_entropy', 'vcg_cdegree_mean', 'vcg_cdegree_variance', 'vcg_cdegree_min', 'vcg_cdegree_max', 'vcg_cdegree_entropy'

  - based on the Variable Graph: 'vg_degree_mean', 'vg_degree_variance', 'vg_degree_min', 'vg_degree_max', 'vg_degree_entropy'

  - based on the Clause Graph: 'cg_degree_mean', 'cg_degree_variance', 'cg_degree_min', 'cg_degree_max', 'cg_degree_entropy'

### Gate Features for CNF Instances

Runs the gate feature extraction algorithm described in~\cite{} and extracts features denoting the numbers of gates of specific types and their distribution of the levels of the detected hierarchical gate structure.

## WCNF Context

### Isohash for WCNF Instances

### Base Features for WCNF Instances

## OPB Context

### Base Features for OPB Instances
