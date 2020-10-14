# peterk87/nf-ionampliseq

Read mapping, variant calling and consensus sequence generation workflow for Ion Torrent Ampliseq sequence data of [FMDV], [CSFV], Zika virus ([ZIKV]), Ebola virus ([EBOV]).

[![GitHub Actions CI Status](https://github.com/peterk87/nf-ionampliseq/workflows/nf-core%20CI/badge.svg)](https://github.com/peterk87/nf-ionampliseq/actions)
[![GitHub Actions Linting Status](https://github.com/peterk87/nf-ionampliseq/workflows/nf-core%20linting/badge.svg)](https://github.com/peterk87/nf-ionampliseq/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/ionampliseq.svg)](https://hub.docker.com/r/nfcore/ionampliseq)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23ionampliseq-4A154B?logo=slack)](https://nfcore.slack.com/channels/ionampliseq)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Workflow Overview

This workflow includes several built-in analysis packages for Ion Torrent AmpliSeq sequence data of [CSFV], [FMDV], [ZIKV] and [EBOV]. Users can also specify their own analysis packages, however, these files must be compatible with the Ion Torrent Software Suite including [tmap] and [tvc]

### Input

- `--sample_sheet`
  - CSV file with 2 columns:
    - Column 1: sample name
    - Column 2: path to raw BAM file from Ion Torrent (absolute path recommended)
- `--org`
  - one of `csf`, `fmd`, `zika` or `ebola` for built-in analysis package, otherwise, user will need to specify a reference genome(s) FASTA file and detailed BED file

### Steps

1. Reference genome selection using [Mash] screen for read mapping
2. Read mapping to top reference genome with [tmap]
3. Variant calling with [tvc]
4. Normalization and filtering of variants using [bcftools]
5. Majority consensus sequence generation using []

[CSFV]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
[FMDV]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=12110&lvl=3&lin=f&keep=1&srchmode=1&unlock
[ZIKV]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=64320&lvl=3&lin=f&keep=1&srchmode=1&unlock
[EBOV]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=186536&lvl=3&lin=f&keep=1&srchmode=1&unlock

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run peterk87/nf-ionampliseq -profile test,<docker/singularity/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run peterk87/nf-ionampliseq -profile <docker/singularity/conda/institute> --input '*_R{1,2}.fastq.gz' --genome GRCh37
    ```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The peterk87/nf-ionampliseq pipeline comes with documentation about the pipeline which you can read at [https://peterk87/nf-ionampliseq/docs](https://peterk87/nf-ionampliseq/docs) or find in the [`docs/` directory](docs).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

## Credits

peterk87/nf-ionampliseq was originally written by Peter Kruczkiewicz.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#ionampliseq` channel](https://nfcore.slack.com/channels/ionampliseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  peterk87/nf-ionampliseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

[tmap]: https://github.com/iontorrent/TS/
[tvc]: http://updates.iontorrent.com/tvc_standalone/
[variantCaller]: https://github.com/iontorrent/TS/tree/master/plugin/variantCaller
