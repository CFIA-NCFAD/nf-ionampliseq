# Multiple sequence alignment (MSA) of FMDV genome sequences

- Input FASTA: `data/fmd/FMDV_AmpliSeq_WG0022620160728referenceSequences.fasta`
- Output MSA FASTA: `data/fmd/msa/fmdv.mafft.fa`

## MAFFT MSA

MSA was generated with [MAFFT][] (v7.471) with the following command:

```bash
$ mafft --thread 4 --auto FMDV_AmpliSeq_WG0022620160728referenceSequences.fasta > fmdv.mafft.fa
```

The strategy selected was FFT-NS-i according to the log:

```
dvtditr (nuc) Version 7.471
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
4 thread(s)


Strategy:
 FFT-NS-i (Standard)
 Iterative refinement method (max. 2 iterations)
```

[MAFFT]: https://mafft.cbrc.jp/alignment/software/