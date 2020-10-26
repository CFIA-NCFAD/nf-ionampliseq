#!/usr/bin/env python

from pathlib import Path

import click
import pandas as pd

mash_screen_columns = ['identity', 'shared_hashes', 'multiplicity', 'pvalue', 'query_id', 'query_comment']

table_front_matter = """
# plot_type: 'table'
# section_name: 'Mash Screen Top Matches'
# section_href: 'https://www.ncbi.nlm.nih.gov/pubmed/31690338'
# description: 'The top matching reference genome for each sample was determined by Mash screen analysis where raw sample read sequences were screened against a Mash sketch (MinHash hash) database of reference genome sequences.'
# pconfig:
#     namespace: 'Mash Screen Results'
# headers:
#     sample:
#         title: 'Sample'
#         description: 'Sample'
#     query_id:
#         title: 'Top Reference'
#         description: 'Top query reference genome based on Mash screen identity'
#     identity:
#         title: 'Identity'
#         description: 'Mash screen identity'
#         format: '{:.1%}'
#     shared_hashes:
#         title: 'Shared Hashes'
#         description: 'Number of shared hashes (or kmers) between query reference genome and sample'
#     multiplicity:
#         title: 'Median Multiplicity'
#         description: 'Mash screen median multiplicity or median number of times shared hashes appear in sample reads'
#     pvalue:
#         title: 'p-value'
#         description: 'Mash screen p-value'
"""


@click.command()
@click.argument('input_directory', type=click.Path(exists=True))
@click.argument('mqc_summary_output', type=click.Path(exists=False))
def main(input_directory,
         mqc_summary_output):
    """Summarize Mash screen results for MultiQC"""
    dir_path = Path(input_directory)
    print(f'input_directory="{dir_path}"')
    dfs = []
    for tsv_path in dir_path.glob('*-mash_screen.tsv'):
        df = pd.read_table(tsv_path, names=mash_screen_columns)
        df = df.sort_values('identity', ascending=False).head(1)
        df['sample'] = tsv_path.stem.replace('-mash_screen', '')
        dfs += [df]
    df_top_results = pd.concat(dfs)
    df_top_results.sort_values('sample', inplace=True)
    df_top_results.set_index('sample', inplace=True)
    with open(mqc_summary_output, 'w') as fout:
        fout.write(table_front_matter)
        df_top_results.to_csv(fout, sep='\t')


if __name__ == '__main__':
    main()
