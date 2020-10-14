#!/usr/bin/env python

import json
from collections import OrderedDict
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
    conseq = str(reccon.seq).upper()
    con_n_chars = conseq.count('N')
    refseq = str(recref.seq).upper()
    ref_n_chars = refseq.count('N')
    aln = edlib.align(conseq, refseq, task='path')
    nice_aln = edlib.getNiceAlignment(aln, conseq, refseq)
    aln['sample'] = sample
    aln['reference'] = recref.id
    aln['consensus_n_chars'] = con_n_chars
    aln['reference_n_chars'] = ref_n_chars
    matches = nice_aln['matched_aligned'].count('|')
    aln['matches'] = matches
    aln['mismatches'] = len(nice_aln['matched_aligned']) - matches
    aln['identity'] = matches / len(conseq)
    if output_aln:
        with open(output_aln, 'w') as fh:
            write_header(fh, 'Input Attributes')
            attrs = OrderedDict()
            attrs['Sample name'] = sample
            attrs['Consensus file'] = consensus
            attrs['Reference file'] = reference
            attrs['Consensus seq name'] = reccon.description
            attrs['Consensus seq length'] = len(conseq)
            attrs['Consensus seq Ns'] = con_n_chars
            attrs['Reference seq name'] = recref.description
            attrs['Reference seq length'] = len(refseq)
            attrs['Reference seq Ns'] = ref_n_chars
            max_attr_key_len = max(len(x) for x in attrs.keys())
            for k, v in attrs.items():
                fh.write(f'{k}:{" " * (max_attr_key_len - len(k))} {v}\n')

            write_header(fh, 'Edlib Alignment')
            aln_keys = ['editDistance', 'alphabetLength', 'locations', 'cigar']
            max_aln_key_len = max(len(x) for x in aln_keys)
            for x in aln_keys:
                fh.write(f'{x}:{" " * (max_aln_key_len - len(x))} {aln[x]}\n')
            write_nice_alignment(nice_aln, conseq, fh, recref, refseq, sample)

    if output_json:
        with open(output_json, 'w') as fh_json:
            json.dump(aln, fh_json)


def write_nice_alignment(nice_aln, conseq, fh, recref, refseq, sample):
    write_header(fh, f'Full alignment - {sample} (Query) VS {recref.id} (Target)')
    fh.write('"M" for match where "." denotes a mismatch and "|" denotes a match')
    aln_len = len(nice_aln['query_aligned'])
    line_width = 60
    for i in range(0, aln_len, line_width):
        for key in ['query_aligned', 'matched_aligned', 'target_aligned']:
            fh.write(f'{key[0].upper()} {i + 1: 5} {nice_aln[key][i:i + line_width + 1]}\n')
        fh.write(f'\n')


def write_header(fh, header):
    fh.write('=' * 80 + '\n' + header + '\n' + '=' * 80 + '\n')


if __name__ == '__main__':
    main()
