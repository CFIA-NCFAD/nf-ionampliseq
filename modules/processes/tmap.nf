process TMAP {
  tag "$sample"

  label 'process_high'

  container 'ghcr.io/cfia-ncfad/nf-ionampliseq:2.0.0'

  input:
  tuple val(sample), path(bam), path(ref_fasta)

  output:
  tuple val(sample), path('*-tmap.{bam,bam.bai}'), path("*-ref.fasta"), emit: bam
  path("versions.yml"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def fasta = "${sample}-ref.fasta"
  def samtools_view_flag = params.output_unmapped_reads ? "" : "-F 4"
  """
  ln -s ${ref_fasta} ${fasta}
  cp ${ref_fasta} ref.fasta
  tmap index -f ref.fasta
  tmap mapall \\
    -f ref.fasta \\
    -r ${bam} \\
    -i bam \\
    -n ${task.cpus} \\
    ${args} \\
    | samtools sort -@${task.cpus} \\
    | samtools view -b ${samtools_view_flag} \\
    > ${sample}-tmap.bam
  samtools index ${sample}-tmap.bam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    tmap: \$(tmap --version 2>&1 | grep 'Version: ' | sed 's/.*Version: //; s/ .*//' | sed 's/\\x1b\\[[0-9;]*[mGK]//g')
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
