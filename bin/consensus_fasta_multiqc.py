#!/usr/bin/env python

import base64
from pathlib import Path

import click


header = '''
<!--
id: 'consensus-fasta'
section_name: 'Consensus sequence FASTA files'
description: 'This section contains the consensus sequences FASTA files for each sample.'
-->

<ul>
'''

footer = '''
</ul>
'''

def to_html_li(fasta: str) -> str:
    fasta_path = Path(fasta)
    b64 = base64.encodebytes(fasta_path.read_bytes())
    return f'''
    <li>
    <a
      download="{fasta_path.name}"
      href="data:text/plain;base64,{b64.decode()}">
        Download {fasta_path.name}
    </a>
    </li>
    '''


@click.command()
@click.argument('input_fasta', type=click.Path(exists=True), nargs=-1)
@click.argument('output_html', type=click.Path(exists=False), nargs=1)
def main(input_fasta, output_html):
    """Embed FASTA file content into HTML for inclusion in MultiQC report"""
    with open(output_html, 'w') as fout:
        fout.write(header)
        for fasta in input_fasta:
            fout.write(to_html_li(fasta))
        fout.write(footer)

if __name__ == '__main__':
    main()
