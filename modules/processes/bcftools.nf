process BCFTOOLS_VCF_NORM {
  tag "$sample"
  publishDir "${params.outdir}/variants/vcf",
             pattern: "*.norm.vcf",
             mode: params.publish_dir_mode

  input:
  tuple val(sample), path(bam), path(vcf), path(ref_fasta)

  output:
  tuple val(sample), path(bam), path('*.norm.vcf'), path(ref_fasta)

  script:
  """
  bcftools norm -Ov -m- -f $ref_fasta $vcf > ${sample}.norm.vcf 
  """
}

process BCFTOOLS_VCF_FILTER {
  tag "$sample"
  publishDir "${params.outdir}/variants/vcf",
             pattern: "*.filt.vcf",
             mode: params.publish_dir_mode

  input:
  tuple val(sample), path(bam), path(vcf), path(ref_fasta)

  output:
  tuple val(sample), path(bam), path('*.filt.vcf'), path(ref_fasta)

  script:
  """
  bcftools filter -e 'INFO/AF<0.5' -Ov $vcf > ${sample}.norm.filt.vcf
  """
}

process BCFTOOLS_STATS {
  tag "$sample"
  publishDir "${params.outdir}/variants/bcftools",
             pattern: "*.bcftools_stats.txt",
             mode: params.publish_dir_mode

  input:
  tuple val(sample), path(bam), path(vcf), path(ref_fasta)

  output:
  path("*.bcftools_stats.txt")

  script:
  """
  bcftools stats -F $ref_fasta $vcf > ${sample}.bcftools_stats.txt
  """
}
