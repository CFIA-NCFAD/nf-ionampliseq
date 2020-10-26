process SAMTOOLS_STATS {
  tag "$sample"
  publishDir "${params.outdir}/samtools/$sample",
             mode: params.publish_dir_mode
  
  input:
  tuple val(sample), path(bam), path(ref_fasta)

  output:
  path('*.{flagstat,idxstats,stats}')

  script:
  """
  samtools flagstat *.bam > ${sample}.flagstat
  samtools idxstats *.bam > ${sample}.idxstats
  samtools stats *.bam > ${sample}.stats
  """
}

process SAMTOOLS_DEPTH {
  tag "$sample"
  publishDir "${params.outdir}/samtools/depth",
             pattern: "*-depths.tsv",
             mode: params.publish_dir_mode

  input:
  tuple val(sample), path(bam), path(ref_fasta)

  output:
  tuple val(sample), path('*-depths.tsv'), path(ref_fasta)

  script:
  """
  samtools depth -aa -d 0 ${bam[0]} > ${sample}-depths.tsv
  """
}
