# peterk87/nf-ionampliseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [2020-10-25]

WIP: First release of peterk87/nf-ionampliseq workflow.

### `Added`

- Edlib pairwise global alignment between consensus and reference sequences. Edlib alignment stats included in table in MultiQC report.
- Updated MultiQC report to by default include more useful information.
- Ion Torrent test dataset at [github.com/peterk87/nf-test-datasets](https://github.com/peterk87/nf-test-datasets/) for use with `-profile test`
- Output documentation

### `Fixed`

- CI and linting Github Actions workflows.
- Documentation

### `Dependencies`

Python:

- [Edlib](https://github.com/Martinsos/edlib) for pairwise alignment and edit distance calculation between strings.
- [pysam](https://pysam.readthedocs.io/en/latest/) for reading BAM file header information.
- [odfpy](https://github.com/eea/odfpy) for reading LibreOffice spreadsheet ODS files into Pandas DataFrames.

## v1.0dev - [date]

Initial release of peterk87/nf-ionampliseq, created with the [nf-core](https://nf-co.re/) template.
