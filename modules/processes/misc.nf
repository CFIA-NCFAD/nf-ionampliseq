/*
* Utility processes
*/

process SAMPLE_INFO_FROM_BAM {
  publishDir "${params.outdir}/bam_sample_info",
             pattern: "*.tsv",
             mode: 'copy'
  input:
  path(bam)

  output:
  tuple path(bam),
        path('sample_name.txt'),
        path('ampliseq_panel.txt'), emit: sample_info
  path '*.tsv'

  script:
  """
  parse_bam_sample_info.py \\
    -i $bam \\
    -o sample_name.txt \\
    -p ampliseq_panel.txt \\
    --write-sample-info
  """
}

process CHECK_SAMPLE_SHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path samplesheet

    output:
    path "samplesheet_reformat.csv" 

    script:
    """
    check_sample_sheet.py $samplesheet samplesheet_reformat.csv
    """
}

process BAM_TO_FASTQ {
  tag "$sample"
  publishDir "${params.outdir}/reads/fastq",
             pattern: "*.fastq.gz",
             mode: 'copy'

  input:
  tuple val(sample), path(bam)

  output:
  tuple val(sample), path("*.fastq.gz")

  script:
  """
  samtools fastq $bam | pigz -c - > ${sample}.fastq.gz
  """
}

process FILTER_BED_FILE {
  tag "$sample"

  input:
  tuple val(sample), path(top_ref), path(bed_file)

  output:
  tuple val(sample), path("*-filtered.bed")

  script:
  """
  REF=\$(head -n1 $top_ref | awk '{print \$5}')
  head -n1 $bed_file > ${sample}-filtered.bed
  grep "^\$REF\\b" $bed_file >> ${sample}-filtered.bed
  """
}
