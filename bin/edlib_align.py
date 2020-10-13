#!/usr/bin/env python

import json
from pathlib import Path


import click
import edlib
import pandas
from Bio import SeqIO


@click.command()
@click.option('--sample', required=True)
@click.option('--consensus', required=True, type=click.Path(exists=True))
@click.option('--reference', required=True, type=click.Path(exists=True))
@click.option('--output-aln', type=click.Path(exists=False))
@click.option('--output-json', type=click.Path(exists=False))
def main(sample,
         consensus,
         reference,
         output_aln,
         output_json):
    """Perform edlib global pairwise alignment of consensus sequence against reference sequence"""
    
    reccon = list(SeqIO.parse(consensus, 'fasta'))[0]
    recref = list(SeqIO.parse(reference, 'fasta'))[0]
    conseq = str(reccon.seq)
    refseq = str(recref.seq)
    aln = edlib.align(conseq, refseq, task='path')
    if output_aln:
        with open(output_aln, 'w') as fh:
            fh.write(f'Sample name: {sample}\n'
                     f'Consensus file: {consensus}\n'
                     f'Reference file: {reference}\n'
                     f'Consensus seq id: {reccon.id}\n'
                     f'Consensus seq length: {len(conseq)}\n'
                     f'Reference seq id: {recref.id}\n'
                     f'Reference seq length: {len(refseq)}\n')
            for x in ['editDistance', 'alphabetLength', 'locations', 'cigar']:
                fh.write(f'edlib alignment {x}: {aln[x]}\n')
            fh.write('='*80 + '\n')
            fh.write(f'Full alignment - {sample} VS {recref.id}\n')
            fh.write('='*80 + '\n')
            nicealn = edlib.getNiceAlignment(aln, conseq, refseq)
            N = len(nicealn['query_aligned'])
            L = 80
            for i in range(0, N, L):
                for x in ['query_aligned', 'matched_aligned', 'target_aligned']:
                    fh.write(f'{i: 5} {nicealn[x][i:i+L+1]}\n')
                fh.write(f'\n')

    if output_json:
        with open(output_json, 'w') as fh_json:
            json.dump(aln, fh_json)

if __name__ == '__main__':
    main()
