#!/usr/bin/env python

import json
from pathlib import Path

import click
import pandas as pd

edlib_columns = [
    ('sample', str),
    ('reference', str),
    ('editDistance', str),
    ('identity', float),
    ('matches', int),
    ('mismatches', int),
    ('locations', object),
    ('consensus_n_chars', int),
    ('reference_n_chars', int),
    ]

table_front_matter = """# plot_type: 'table'
# section_name: 'Edlib Global Pairwise Alignment'
# section_href: 'https://github.com/Martinsos/edlib'
# description: 'Global pairwise alignment of read mapping/variant calling majority consensus sequence against reference sequence using Edlib.'
# pconfig:
#     namespace: 'Edlib'
# headers:
#     sample:
#         title: 'Sample'
#         description: 'Sample name'
#         format: '{}'
#     reference:
#         title: 'Reference'
#         description: 'Reference sequence ID'
#         scale: False
#         format: '{}'
#     editDistance:
#         title: 'Edit Distance'
#         description: 'Edlib edit (Levenshtein) distance'
#         format: '{:.0f}'
#     identity:
#         title: 'Identity'
#         description: 'Mash screen identity'
#         format: '{:.1%}'
#     matches:
#         title: 'Matches'
#         description: 'Edlib alignment matching positions'
#         format: '{:.0f}'
#     mismatches:
#         title: 'Mismatches'
#         description: 'Edlib alignment mismatching positions'
#         format: '{:.0f}'
#     locations:
#         title: 'Locations'
#         description: 'Edlib alignment locations'
#         scale: False
#         format: '{}'
#     consensus_n_chars:
#         title: 'Ns Sample'
#         description: 'Number of Ns in sample consensus sequence'
#         format: '{:.0f}'
#     reference_n_chars:
#         title: 'Ns Ref'
#         description: 'Number of Ns in reference sequence'
#         format: '{:.0f}'
"""


@click.command()
@click.argument('input_directory', type=click.Path(exists=True))
@click.argument('mqc_summary_output', type=click.Path(exists=False))
def main(input_directory,
         mqc_summary_output):
    """Summarize Edlib pairwise alignment results for MultiQC"""
    dir_path = Path(input_directory)
    res = []
    for p in dir_path.glob('*.json'):
        with open(p) as f:
            res.append(json.load(f))
    df = pd.DataFrame(res)[[x for x, y in edlib_columns]]
    df = df.astype({x:y for x, y in edlib_columns})
    print(df)
    df.sort_values('sample', inplace=True)
    df.set_index('sample', inplace=True)
    with open(mqc_summary_output, 'w') as fout:
        fout.write(table_front_matter)
        df.to_csv(fout, sep='\t', quoting=None)


if __name__ == '__main__':
    main()
