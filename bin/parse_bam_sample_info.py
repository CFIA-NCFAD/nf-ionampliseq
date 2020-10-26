#!/usr/bin/env python

import re
import logging
from collections import OrderedDict

import click
import pysam
import pandas as pd


logger = logging.getLogger(__name__)

regex_fasta_filename = re.compile(r'.+\/([^\/]+)\.fasta\b.*')


@click.command()
@click.option('-i', '--input-bam', type=click.Path(exists=True), required=True)
@click.option('-o', '--output-sample-name', type=click.Path(exists=False), required=True)
@click.option('-p', '--output-ampliseq-panel', type=click.Path(exists=False), required=True)
@click.option('--write-sample-info', is_flag=True, help='Write sample info tab-delimited table')
def main(input_bam,
         output_sample_name,
         output_ampliseq_panel,
         write_sample_info):
    """
    Parse sample info from an Ion Torrent tmap BAM file

    The main info parsed from the BAM file include the sample name and AmpliSeq panel used.
    """
    bam = pysam.AlignmentFile(input_bam)
    sample_info = None
    try:
        rg = bam.header['RG'][0]
        sample_name = rg['SM']
        
        with open(output_sample_name, 'w') as fout:
            fout.write(sample_name)
        pg = bam.header.get('PG', None)
        ref_fasta = None
        panel = ''
        # get filename reference fasta used by tmap
        if pg:
            for x in pg:
                if x['ID'] == 'tmap':
                    cl = x['CL']
            ref_fasta = regex_fasta_filename.sub(r'\1', cl)
            
            if ref_fasta:
                ref_fasta_lc = ref_fasta.lower()
                if ref_fasta_lc.startswith('fmd'):
                    panel = 'fmd'
                elif ref_fasta_lc.startswith('csf'):
                    panel = 'csf'
                else:
                    logger.warning(f'Could not determine AmpliSeq panel of BAM "{input_bam}", sample "{sample_name}". Reference fasta filename was "{ref_fasta}".')
        if output_ampliseq_panel:
            with open(output_ampliseq_panel, 'w') as fout:
                fout.write(panel)
        if write_sample_info:
            sample_info = OrderedDict(
                Sample=sample_name,
                BAM=input_bam,
                AmpliSeq_Panel=panel,
                Reference_FASTA=ref_fasta,
                Platform=rg.get('PL', None),
                Platform_Unit=rg.get('PU', None),
                Run_Date=rg.get('DT', None),
                ID=rg.get('ID', None),
                Seq_Centre=rg.get('CN', None),
                Program_Group=rg.get('PG', None),
                Key_Sequence=rg.get('KS', None),
                Flow_Order=rg.get('FO', None),)
            df = pd.DataFrame([sample_info])
            df.to_csv(f'{sample_name}-bam_sample_info.tsv', index=False, sep='\t')
    except KeyError as e:
        raise KeyError(f'Could not find sample name in BAM header of "{input_bam}". First read group (RG) entry was {rg}')


if __name__ == '__main__':
    main()