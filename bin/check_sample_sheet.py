#!/usr/bin/env python

from pathlib import Path

import click
import pandas as pd


@click.command()
@click.argument('input_sample_sheet', type=click.Path(exists=True))
@click.argument('output_sample_sheet', type=click.Path(exists=False))
def main(input_sample_sheet,
         output_sample_sheet):
    """Check and reformat sample sheet into CSV"""
    input_path = Path(input_sample_sheet)
    print(f'input_sample_sheet="{input_path}"')
    ext = input_path.suffix.lower()
    print(f'ext={ext}')
    try:
        if ext in ['.tsv', '.txt', '.tab']:
            df = pd.read_table(input_path)
        elif ext == '.csv':
            df = pd.read_csv(input_path)
        elif ext in ['.xls', '.xlsx', '.ods']:
            df = pd.read_excel(input_path)
        else:
            raise ValueError(f'Unknown file format for sample sheet "{input_path}"')
    except Exception as ex:
        print(ex)
        raise ex
    print(df)
    assert df.shape[1] == 2, f'Two columns expected in sample sheet, but {df.shape[1]} found!'
    df.columns = ['sample', 'bam_filepath']
    bam_filepaths = []
    for fp in df.bam_filepath:
        path = Path(fp)
        assert path.suffix.lower() == '.bam', \
            f'BAM file path extension is "{path.suffix.lower()}" not ".bam"! ' \
            f'Please check that you have specified the correct file.'
        bam_filepaths.append(str(path.resolve().absolute()))
    df.bam_filepath = bam_filepaths
    df.to_csv(output_sample_sheet, index=False)
    print(f'Wrote reformatted sample sheet CSV to "{output_sample_sheet}"')


if __name__ == '__main__':
    main()
