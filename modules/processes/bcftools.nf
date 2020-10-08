process BCFTOOLS_VCF_NORM {
  tag "$sample"
  publishDir "${params.outdir}/variants/vcf",
             pattern: "*.norm.vcf",
             mode: 'copy'

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
             mode: 'copy'

  input:
  tuple val(sample), path(bam), path(vcf), path(ref_fasta)

  output:
  tuple val(sample), path(bam), path('*.filt.vcf'), path(ref_fasta)

  script:
  """
  bcftools filter -e 'INFO/AF<0.5' -Ov $vcf > ${sample}.norm.filt.vcf
  """
}
