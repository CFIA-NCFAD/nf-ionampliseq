# CFIA-NCFAD/nf-ionampliseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[2.0.1](https://github.com/CFIA-NCFAD/nf-ionampliseq/releases/tag/2.0.1))] - 2025-08-13

This release makes the potential frameshift introducing filter for variant calling results optional and off by default.

### Added

- [config] `filter_frameshift_variants` pipeline parameter (default: `false`) for optionally filtering potential frameshift introducing variants from the TVC variant calling results prior to consensus sequence generation.

### Changes

- [docs] Updated `README.md` and citation information.
- [update] Use TVC Flow Evaluator calculated allele fraction (AF) for computing GT in Bcftools filtered VCF instead of observation based counts.

### Fixed

- [fix] Turned off potential frameshift variant filtering leading to inaccurate consensus sequences.
- [fix] Removed git diff action from `docker.yml` since features are unsupported by GHA.

## [[2.0.0](https://github.com/CFIA-NCFAD/nf-ionampliseq/releases/tag/2.0.0))] - 2025-08-12

This release modernizes nf-ionampliseq to DSL-2 standards with enhanced BLAST analysis capabilities, improved containerization, and better resource management. Key additions include BLAST database integration, enhanced variant calling options, and comprehensive process version tracking.

### Added

- [feat] `BLASTN` process to run `blastn` with consensus sequences against user-specifed DB (`--blast_db`).
- [feat] `BLASTN_COVERAGE` process and `blast_coverage.py` for summarizing BLAST results and generating MultiQC tables.
- [feat] `fill-tags` plugin for `BCFTOOLS_FILTER` and better filtering of variants for consensus sequence construction.
- [feat] `SEQTK_SUBSEQ` process for extracting reference sequences based on Mash screen results.
- [feat] `CONSENSUS_MULTIQC` process for enhanced consensus sequence reporting in MultiQC.
- [feat] Enhanced container support for Docker, Singularity, Apptainer, and Podman.
- [feat] Version tracking for all processes with `versions.yml` files.
- [config] comprehensive DSL2 module configuration in `conf/modules.config`.
- [config] new configuration parameters: `blast_db`, `trim_primers`, `output_unmapped_reads`.
- [config] variant calling options: `minor_allele_fraction`, `major_allele_fraction`, `low_coverage`.
- [config] enhanced container profiles for multiple container engines.
- [config] proper resource management with `check_max` function for memory, CPU, and time limits.
- [config] timestamped execution reports and timeline files in `nextflow.config`.

### Changes

- [docs] Comprehensive documentation overhaul with detailed parameter descriptions and external tool references
- [cleanup] Removed legacy plotting module and environment.yml file for simplified architecture. Use wgscovplot instead for better and interactive plotting of coverage stats and variants.
- [cleanup] Simplified TMAP process configuration using dynamic arguments parameter.
- [update] Repository and container registry references updated to CFIA-NCFAD organization
- [feat] Complete workflow modernization to DSL2 standards with proper process separation.
- [feat] Added new Docker build workflow in `.github/workflows/docker.yml` for improved container publishing and CI separation.
- [feat] Added support for building and using custom Docker images with updated samtools and TMAP/TVC binaries.
- [feat] Enhanced process version tracking: all major processes now emit `versions.yml` with tool versions for reproducibility.
- [feat] Improved error handling and retry strategies for all processes.
- [update] MultiQC process updated to support new BLAST coverage MultiQC table.
- [update] Bump version of workflow and dependencies for 2.0.0 release.
- [update] Repository renamed from `peterk87/nf-ionampliseq` to `CFIA-NCFAD/nf-ionampliseq`.
- [update] Nextflow version requirement updated to `!>=22.10.1`.
- [update] All processes now use proper container definitions with conda and container specifications.
- [update] Process resource labeling standardized (`process_low`, `process_medium`, `process_high`).
- [update] Input/output patterns standardized across all processes with proper emit declarations.
- [update] TMAP and TVC processes enhanced with better parameter handling and version tracking.
- [update] TMAP and TVC processes now use new Docker image `ghcr.io/cfia-ncfad/nf-ionampliseq:2.0.0` with updated samtools and runtime dependencies.
- [update] MASH screen workflow restructured with separate sketching and screening processes.
- [update] FastQC process enhanced with memory optimization and version tracking.
- [update] Samtools processes standardized with consistent version tracking and output handling.
- [update] Mosdepth process optimized with improved output handling and version tracking.
- [update] Edlib processes enhanced with container support and version tracking.
- [update] Sample sheet processing improved with better error handling and validation.
- [docs] Documentation and help text for BLAST coverage analysis improved.
- [docs] Updated documentation to reflect new container build and usage instructions.
- [ci] Updated GitHub Actions CI tests.
- [ci] Moved Docker container build to separate workflow for improved CI/CD.
- [ci] CI now tests with multiple Nextflow versions 22.10.1 and latest stable (25.04.6 currently).

### Fixed

- [fix] Container compatibility issues across different container engines.
- [fix] Process resource allocation and memory management.
- [fix] Version tracking consistency across all processes.
- [fix] Input/output file handling and validation.
- [fix] MultiQC report generation and custom content integration.

### Dependencies

- [deps] Updated to Nextflow `!>=22.10.1`.
- [deps] Enhanced container support with multiple engine options.
- [deps] Standardized conda environment specifications across all processes.
- [deps] Added better version tracking for all bioinformatics tools.

## v1.0.1 - [2020-11-06]

Patch release to make MultiQC reports a bit more useful with consensus sequence FASTA files embedded for download.

### Added

- Consensus sequence FASTA files are embedded within the MultiQC HTML report and can be downloaded from it.

## v1.0.0 - [2020-10-25]

First release of peterk87/nf-ionampliseq workflow.

### Added

- Edlib pairwise global alignment between consensus and reference sequences. Edlib alignment stats included in table in MultiQC report.
- Updated MultiQC report to by default include more useful information.
- Ion Torrent test dataset at [github.com/peterk87/nf-test-datasets](https://github.com/peterk87/nf-test-datasets/) for use with `-profile test`
- Output documentation

### Fixed

- CI and linting Github Actions workflows.
- Documentation

### Dependencies

Python:

- [Edlib](https://github.com/Martinsos/edlib) for pairwise alignment and edit distance calculation between strings.
- [pysam](https://pysam.readthedocs.io/en/latest/) for reading BAM file header information.
- [odfpy](https://github.com/eea/odfpy) for reading LibreOffice spreadsheet ODS files into Pandas DataFrames.

## v1.0dev - [date]

Initial release of peterk87/nf-ionampliseq, created with the [nf-core](https://nf-co.re/) template.
