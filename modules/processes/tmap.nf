process TMAP {
  tag "$sample"

  label 'process_high'

  container 'peterk87/nf-ionampliseq:1.0.0'

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
  ln -s $ref_fasta $fasta
  cp $ref_fasta ref.fasta
  tmap index -f ref.fasta
  tmap mapall \\
    -f ref.fasta \\
    -r $bam \\
    -i bam \\
    -o 2 \\
    -n ${task.cpus} \\
    -u -v --do-realign \\
    --prefix-exclude 5 \\
    -Y -J 25 --end-repair 15 \\
    --do-repeat-clip \\
    --context stage1 map4 \\
    | samtools sort -@${task.cpus} \\
    | samtools view -b $samtools_view_flag \\
    > ${sample}-tmap.bam
  samtools index ${sample}-tmap.bam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    tmap: \$(tmap --version 2>&1 | grep '^Version: ' | sed 's/^Version: //; s/ .*\$//')
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}


// TODO: fix:
// tmap: \$(tmap --version 2>&1 | grep '^version: ' | sed 's/^Version: //; s/ .*\$//')
// use ${args} from modules.config instead of the static args
