# FMDV Phylogenetic Tree Construction

A maximum-likelihood phylogenetic tree was constructed from a multiple sequence alignment (MSA) of 457 FMDV genome sequences using [IQ-TREE][] (v2.0.3).

Best model selected by ModelFinder was GTR+F+R9:

```
Akaike Information Criterion:           GTR+F+R10
Corrected Akaike Information Criterion: GTR+F+R9
Bayesian Information Criterion:         GTR+F+R9
Best-fit model: GTR+F+R9 chosen according to BIC

All model information printed to fmdv.mafft.fa.model.gz
CPU time for ModelFinder: 33412.808 seconds (9h:16m:52s)
Wall-clock time for ModelFinder: 4198.642 seconds (1h:9m:58s)
```

## IQ-TREE Version

```bash
 $ iqtree --version
IQ-TREE multicore version 2.0.3 for Linux 64-bit built Apr 26 2020
Developed by Bui Quang Minh, Nguyen Lam Tung, Olga Chernomor,
Heiko Schmidt, Dominik Schrempf, Michael Woodhams.
```

[IQ-TREE]: http://www.iqtree.org/
