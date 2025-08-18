# CFIA-NCFAD/nf-ionampliseq: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [BAM Sample Info](#bam-sample-info) - Sample info extracted from BAM file headers
* [FASTQ Reads](#fastq-reads) - BAM to FASTQ output
* [FastQC](#fastqc) - Read quality control
* [Mash](#mash) - Top reference genome determination by [Mash][] screen with configurable k-mer size and sketch count
* [TMAP](#tmap) - Read mapping using the Thermo Fisher mapper [tmap] with optional unmapped reads output
* [Samtools](#samtools) - Read mapping stats calculation with [Samtools][]
* [Mosdepth](#mosdepth) - Coverage stats calculated by [Mosdepth][]
* [TVC](#tvc) - Variant calling using the Thermo Fisher variant caller [tvc] with AmpliSeq primer trimming and configurable parameters
* [Bcftools](#bcftools) - Variant filtering for majority consensus sequence generation and variant statistics for MultiQC report
* [Consensus Sequence](#consensus-sequence) - Majority consensus sequence with configurable coverage masking thresholds
* [Edlib Pairwise Alignment](#edlib-pairwise-alignment) - Pairwise global alignment and edit distance between reference and consensus sequences
* [BLAST Analysis](#blast-analysis) - BLAST database search against consensus sequences for taxonomic identification and sequence similarity analysis
* [CSFV Consensus Sequence QC](#csfv-consensus-sequence-qc) - Quality control analysis specific to CSFV consensus sequences
* [MultiQC](#multiqc) - Aggregate report describing results from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## BAM Sample Info

Sample and sequencing run information is extracted from the BAM file headers and with certain pipeline run modes (`--input`), it is used to determine which AmpliSeq panel to use and what the sample names are for each BAM file.

**Output files:**

* `bam_sample_info/`
  * `*.tsv`: Sample info table containing fields such as platform, platform unit, run date, AmpliSeq panel, original reference FASTA, sample name, etc.

## FASTQ Reads

FASTQ format read sequences are extracted from input BAM files.

**Output files:**

* `reads/`
  * `fastq/`
    * `*.fastq.gz`: Gzip compressed FASTQ format read sequences extracted from input BAM files for each sample. *Useful for troubleshooting and ad hoc analysis.*

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output files:**

* `fastqc/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.
* `fastqc/zips/`
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots displayed in the MultiQC report shows *untrimmed* reads. They may contain adapter sequence and potentially regions with low quality.

## Mash

[Mash][] [screen](https://doi.org/10.1186/s13059-019-1841-x) is used for determining the top reference genome for each sample. A Mash [MinHash](https://en.wikipedia.org/wiki/MinHash) sketch database of all reference genome sequences is "screened" against a sample's reads to determine "containment" of reference genome sequence (i.e. how many kmers are shared between the references and the reads). The top reference genome by Mash screen identity is used for read mapping, variant calling and majority consensus sequence generation.

**Output files:**

* `mash_screen/`
  * `*-mash_screen-top_ref.fasta`: Top reference genome sequence in FASTA format.
  * `*-mash_screen.tsv`: Raw sorted Mash screen results of sample's reads against all reference genomes showing: identity, shared hashes, median multiplicity, p-value, reference genome ID, and comment (if applicable).

> **NB:** A Mash screen section and results summary table can be found in the MultiQC report.

## TMAP

The Thermo Fisher Scientific read mapper ([TMAP][]) is used for read mapping of the Ion Torrent reads against the top reference genome. The run parameters for [TMAP] are similar to those of the [variantCaller] plugin of the [Torrent Suite] Ion Torrent analysis platform.

**Output files:**

* `tmap/`
  * `*-ref.fasta`: Top reference genome sequence FASTA.
  * `*-tmap.bam`: BAM file with read alignment of sample reads against top reference.
  * `*-tmap.bam.bai`: BAM file index (`samtools index` output).

> **NB:** Due to issues with using FASTQ files as input for TMAP, only the original Ion Torrent BAM files are supported with this workflow.

## Samtools

[Samtools][] is used to compute read mapping statistics and metrics including coverage depth at each position of the reference genome.

**Output files:**

* `samtools/`
  * `${sample}/`
    * `*.flagstat`: The number of alignments for each FLAG type.
    * `*.idxstats`: Alignment summary statistics.
    * `*.stats`: Comprehensive read alignment statistics.

> **NB:** Stats and metrics computed by Samtools are used to generate multiple fields in the MultiQC general statistics table and to generate multiple plots under the section "VARIANTS: Samtools (tmap raw)"

## Mosdepth

[Mosdepth][] is used to compute coverage depth and breadth statistics. Mosdepth results for each sample are visualized and summarized in the MultiQC report.

**Output files:**

* `mosdepth/`
  * `*.mosdepth.global.dist.txt`: contains, a cumulative distribution indicating the proportion of total bases  that were covered for at least a given coverage value. It does this for each chromosome, and for the whole genome.
  * `*.mosdepth.region.dist.txt`: contains, a cumulative distribution indicating the proportion of 200 base window regions that were covered for at least a given coverage value. It does this for each chromosome, and for the whole genome.
  * `*.mosdepth.summary.txt`: summary read alignment depth statistics (reference sequence length, total bases, mean/min/max).
  * `*.per-base.bed.gz`: 4 column BED file with per base coverage. Columns in order contain: reference ID; 0-base start index; 0-base end index; coverage depth.
  * `*.regions.bed.gz`: 4 column BED file with 200 base window regions mean coverage. Columns in order contain: reference ID; 0-base start index; 0-base end index; coverage depth.
  * `*.gz.csi`: Table index files for rapid random read access.

## TVC

[TVC][] performs variant calling from the [tmap][] read mapping output. It will also trim AmpliSeq primers given the coordinates of the amplicons specified in each AmpliSeq panel's detailed BED file. TVC will also account for common Ion Torrent sequencing specific error motifs.

**Output files:**

* `tvc/`
  * `${sample}/`
    * `tvc-outdir`: TVC output directory
      * `*-small_variants.vcf`: SNPs, MNPs and small indels detected by TVC in TMAP read alignment.
      * `*-small_variants_filtered.vcf`: Small variants that were filtered out (TVC help unclear what this file is supposed to contain).
      * `tvc_metrics.json`: JSON file with basic transition/transversion (Ts/Tv) stats.
      * `indel_assembly.vcf`: *De novo* assembly of large indels detected by TVC with SPAdes. This VCF file should typically contain no variants.
      * `black_listed.vcf`: Variants detected that are "black listed" (i.e. variants that one wishes to ignore). This file should be empty in most cases.
    * `*-tvc-postprocessed.bam`: Post-processed TVC output BAM file for debugging purposes.
    * `*-tvc-postprocessed.bam.bai`: Post-processed TVC output BAM index file.
    * `fasta/`
      * `*-ref.fasta`: Top reference genome sequence FASTA selected by Mash screen.

## Bcftools

Variants detected by [TVC] are normalized and filtered by [Bcftools][] to determine the high-confidence minor alleles (allele fraction (AF) >= 0.25 by default) and high-confidence major alleles (AF >= 0.75 by default). Supplemental information is added to the VCF using the Bcftools [fill-tags plugin](https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html). Bcftools is also used to report the number of SNPs, MNPs and indels detected as well as transitions/transversions (Ts/Tv). Bcftools statistics are reported in the MultiQC report general stats table.

**Output files:**

* `variants/`
  * `bcftools_stats/`
    * `${sample}.bcftools_stats.txt`: Variant statistics file used by MultiQC for reporting and visualization.
  * `${sample}/`
    * `${sample}.bcftools_filt.vcf`: Bcftools normalized, supplemented with additional allele info (Bcftools fill-tags plugin) and filtered VCF. This is the VCF that is used for consensus sequence generation with Bcftools. Multiallelic sites are split into individual records for easier filtering and consensus sequence generation.
    * `${sample}-ref.fasta`: Reference sequence used for variant calling.

## BLAST Analysis

[BLAST] is used to search consensus sequences against user-specified databases for virus strain identification and sequence similarity analysis. This step is optional and requires the `--blast_db` parameter to be specified (e.g. `nextflow run CFIA-NCFAD/nf-ionampliseq ... --blast_db /path/to/blast_db/core_nt`).

**Output files:**

* `blast/`
  * `*.blastn.tsv`: BLASTN output in tabular format showing query coverage, identity, and alignment details.

> **NB:** Summarized BLAST analysis results are integrated into the MultiQC report, providing comprehensive coverage statistics and virus strain identification information.

## Consensus Sequence

A majority consensus sequence is constructed with [Bcftools][] `consensus` from a coverage depth masked reference sequence (low/no coverage reference sequence positions are replaced with `N`) and variants normalized and filtered by [Bcftools][] according to the minor and major allele fraction thresholds (workflow params `minor_allele_fraction = 0.25` and `major_allele_fraction = 0.75` by default). Positions with allele frequencies between the minor and major allele fraction threshold will appear as ambiguous nucleotides.

**Output files:**

* `consensus/`
  * `*.consensus.fasta`: Depth-masked consensus sequence. Low/no coverage positions masked with `N` (by default).
  * `vcf/`
    * `*.merged.vcf.gz`: Re-merged split variants used for consensus sequence construction.

## CSFV Consensus Sequence QC

The CSFV consensus sequence quality control process analyzes consensus FASTA files to determine the quality of consensus sequences. This step provides comprehensive quality metrics including N character distribution, GC content, reading frame analysis, and assembly quality assessment. Results are integrated into the MultiQC report for comprehensive quality reporting.

**Output files:**

* `qc_csfv_fasta/`
  * `csfv_quality_summary.csv`: CSV file containing detailed quality metrics for each consensus sequence including coverage percentage, N character distribution, GC content, frameshift evidence, and polyprotein analysis.
  * `csfv_quality_summary_mqc.txt`: MultiQC-compatible table format for integration into the final MultiQC report.

> **NB:** CSFV QC results are automatically integrated into the MultiQC report, providing comprehensive quality assessment alongside other pipeline metrics.

## Edlib Pairwise Alignment

[Edlib][] is used to perform global pairwise alignment between the majority consensus sequence and the top reference sequence for each sample. Edlib reports the edit or Levenshtein distance between the two sequences (minimum number of changes to change one sequence into the other). Edlib results are reported in the MultiQC report. A naive base by base pairwise similarity identity is also calculated for the MultiQC report.

**Output files:**

* `consensus/`
  * `edlib/`
    * `*.edlib.txt`: Statistics and human-readable Edlib pairwise global alignment between reference and consensus sequences. This file may be useful for troubleshooting.
    * `*.json`: Edlib JSON output.

## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`  
  * `multiqc_report.html`: a standalone MultiQC HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats (PNG, PDF, SVG).

> **NB:** All consensus sequence FASTA files will be embedded within the MultiQC HTML report and can be downloaded from it.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.

<!-- External links and references -->

[Mosdepth]: https://github.com/brentp/mosdepth
[BLAST]: https://blast.ncbi.nlm.nih.gov/
[Samtools]: https://www.htslib.org/
[TMAP]: https://github.com/iontorrent/TS/
[TVC]: http://updates.iontorrent.com/tvc_standalone/
[variantCaller]: https://github.com/iontorrent/TS/tree/master/plugin/variantCaller
[Torrent Suite]: https://github.com/iontorrent/TS
[Mash]: https://mash.readthedocs.io/en/latest/
[Bcftools]: https://samtools.github.io/bcftools/bcftools.html
[Edlib]: https://github.com/Martinsos/edlib
