# CFIA-NCFAD/nf-ionampliseq: Usage

## Introduction

This pipeline performs read mapping and variant calling with Thermo Fisher developed open-source tools, [TMAP] and [TVC] to produce more accurate read mapping and variant calling results from Ion Torrent AmpliSeq sequence data. The pipeline also generates a consensus sequence and comprehensive QC stats and results report using MultiQC.

This workflow currently includes several built-in analysis packages for Ion Torrent AmpliSeq sequence data of [CSFV] and [FMDV]. Users can also specify their own AmpliSeq panels, however, these files (reference sequences FASTA and detailed BED file) must be compatible with the Ion Torrent Software Suite including [TMAP] and [TVC].

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run CFIA-NCFAD/nf-ionampliseq --input '/path/to/ion-torrent/*.bam' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Help and usage info

The help and usage information for this pipeline can be displayed in your terminal with:

```bash
nextflow run CFIA-NCFAD/nf-ionampliseq --help
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull CFIA-NCFAD/nf-ionampliseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [CFIA-NCFAD/nf-ionampliseq releases page](https://github.com/CFIA-NCFAD/nf-ionampliseq/releases) and find the latest version number - numeric only (eg. `2.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Pipeline parameters

### Input/output options

Define where the pipeline should find input data and save output data.

#### `--input`

- **Required**
- Type: string

Input BAM files from Ion Torrent Torrent Suite tmap.

Use this to specify the location of your input BAM files. For example:

```bash
--input 'path/to/data/*.bam'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. The BAM files must be produced by the Ion Torrent Torrent Suite software, specifically `tmap`

#### `--rundir`

- **Optional**
- Type: string

Ion Torrent sequencing run directory.

Use this to specify the exported Ion Torrent sequencing run you wish to analyze. For example:

```bash
--rundir path/to/iontorrent/run
```

Please note the following requirements:

1. The directory specified must contain BAM files matching the pattern `IonCode_*_rawlib.bam` and a `ion_params_00.json` file
2. The BAM files must be produced by the Ion Torrent Torrent Suite software, specifically `tmap`

#### `--sample_sheet`

- **Optional**
- Type: string

Sample sheet table containing two fields: sample name, BAM file path.

Use this option to specify the Ion Torrent sequencing samples that you wish to analyze and if you wish to rename samples rather than rely on the names within the BAM files. Using a samplesheet is necessary if you wish to merge the reads of the same sample with different BAM sample names.

Please note the following requirements:

1. The sample sheet table file may be one of the following formats: tab-delimited (TSV), CSV, ODS or XLSX
2. The sample sheet table must contain a header row and 2 fields: the first field must contain the sample name; and the second field must contain the BAM file path or URL
3. The BAM files must be produced by the Ion Torrent Torrent Suite software, specifically `tmap`

#### `--panel`

- **Optional**
- Type: string
- Default: null

AmpliSeq panel ("fmd" or "csf").

Specify a built-in AmpliSeq panel (reference genome sequences FASTA and detailed BED file of AmpliSeq amplicon coordinates) to perform analysis with.

#### `--ref_fasta`

- **Optional**
- Type: string
- Default: null

Reference genome sequences FASTA for AmpliSeq panel.

#### `--bed_file`

- **Optional**
- Type: string
- Default: null

Detailed BED file containing AmpliSeq amplicon coordinates and information.

#### `--output_unmapped_reads`

- **Optional**
- Type: boolean
- Default: false

Whether or not to include unmapped reads in TMAP BAM output.

#### `--blast_db`

- **Optional**
- Type: string
- Default: null

BLAST database prefix for sequence similarity searches (e.g. `nt` or `core_nt`).

Path to a BLAST database (e.g. `nt` or `core_nt`).

#### `--outdir`

- **Optional**
- Type: string
- Default: `./results`

The output directory where the results will be saved.

### Mash Screen Parameters

Mash screen parameters for determining top reference genome from MinHash kmer containment of reference hashes in sample read sequences.

#### `--mash_k`

- **Optional**
- Type: integer
- Default: 19

Reference genome sequence Mash sketch k-mer size.

#### `--mash_s`

- **Optional**
- Type: integer
- Default: 10000

Number of Mash sketches to create for each reference genome sequence.

### TVC Parameters

Thermo Fisher variant caller (TVC) parameters.

#### `--tvc_error_motifs_dir`

- **Optional**
- Type: string
- Default: `$baseDir/data/tvc-sse`

Directory with Ion Torrent TVC error motifs.

#### `--tvc_read_limit`

- **Optional**
- Type: integer
- Default: 4000000

TVC read limit.

#### `--tvc_downsample_to_coverage`

- **Optional**
- Type: integer
- Default: 8000

TVC downsample to at most X coverage.

#### `--tvc_min_mapping_qv`

- **Optional**
- Type: integer
- Default: 0

TVC min mapping quality value.

#### `--tvc_read_snp_limit`

- **Optional**
- Type: integer
- Default: 20

TVC: do not use reads with number of SNPs above this limit.

#### `--trim_primers`

- **Optional**
- Type: boolean
- Default: true

Enable TVC AmpliSeq primer trimming.

### Variant Calling Options

Various options for the variant calling branch of the pipeline.

#### `--minor_allele_fraction`

- **Optional**
- Type: number
- Default: 0.25

Minor variant allele frequency/fraction.

#### `--major_allele_fraction`

- **Optional**
- Type: number
- Default: 0.75

Major variant allele frequency/fraction. Only major variant alleles are used for generating a consensus sequence.

#### `--filter_frameshift_variants`

- Type: boolean
- Default: false

Filter any variants from TVC VCF that may be frameshift mutations (e.g. REF='A', ALT='AG'). This may be too aggressive for Ion Torrent, but usually necessary for ONT data. Inframe indels are preserved. This option may affect the consensus sequence.

#### `--low_coverage`

- **Optional**
- Type: integer
- Default: 1

Low coverage depth threshold. Consensus sequence positions with less than this coverage depth will be masked with the `low_cov_char` character.

#### `--no_coverage`

- **Optional**
- Type: integer
- Default: 0

No coverage depth threshold. Positions at or below this coverage will be masked with the `no_cov_char` character.

#### `--low_cov_char`

- **Optional**
- Type: string
- Default: 'N'

Mask low coverage positions in reference with this character.

#### `--no_cov_char`

- **Optional**
- Type: string
- Default: 'N'

Mask no coverage positions in reference with this character.

### Slurm Scheduler Options

#### `--slurm_queue`

- **Optional**
- Type: string

Slurm queue/partition to submit pipeline jobs to.

#### `--slurm_queue_size`

- **Optional**
- Type: integer
- Default: 100

Slurm queue size. Max number of jobs to queue at once.

### Max Job Request Options

Set the top limit for requested resources for any single job.

#### `--max_cpus`

- **Optional**
- Type: integer
- Default: 4

Maximum number of CPUs that can be requested for any single job.

> **NOTE:** Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`

#### `--max_memory`

- **Optional**
- Type: string
- Default: `32.GB`

Maximum amount of memory that can be requested for any single job.

> **NOTE:** Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`

#### `--max_time`

- **Optional**
- Type: string
- Default: `240.h`

Maximum amount of time that can be requested for any single job.

> **NOTE:** Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`

### Generic Options

Less common options for the pipeline, typically set in a config file.

#### `--help`

- **Optional**
- Type: boolean

Display help text.

#### `--publish_dir_mode`

- **Optional**
- Type: string
- Default: `copy`

Method used to save pipeline results to output directory.

> **NOTE:** The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.

#### `--monochrome_logs`

- **Optional**
- Type: boolean

Do not use coloured log outputs.

> **NOTE:** Set to disable colourful command line output and live life in monochrome.

#### `--tracedir`

- **Optional**
- Type: string
- Default: `${params.outdir}/pipeline_info`

Directory to keep pipeline Nextflow logs and reports.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`CFIA-NCFAD/nf-ionampliseq`](https://hub.docker.com/r/CFIA-NCFAD/nf-ionampliseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`CFIA-NCFAD/nf-ionampliseq`](https://hub.docker.com/r/CFIA-NCFAD/nf-ionampliseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `TVC` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: TVC {
    memory = 32.GB
  }
}
```

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<!-- External links and references -->

[CSFV]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=3052625&lvl=3&keep=1&srchmode=1&unlock
[FMDV]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=12110&lvl=3&lin=f&keep=1&srchmode=1&unlock
[TMAP]: https://github.com/iontorrent/TS/
[TVC]: http://updates.iontorrent.com/tvc_standalone/
