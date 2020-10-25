process TMAP {
  tag "$sample"
  publishDir "${params.outdir}/tmap",
             pattern: "*-tmap.bam*",
             mode: params.publish_dir_mode
  publishDir "${params.outdir}/tmap",
             pattern: "*.fasta",
             saveAs: { "$sample-ref.fasta" },
             mode: params.publish_dir_mode

  input:
  tuple val(sample), path(bam), path(ref_fasta)

  output:
  tuple val(sample), path('*-tmap.{bam,bam.bai}'), path(ref_fasta)

  script:
  """
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
    | samtools sort -O BAM -o ${sample}-tmap.bam
  samtools index ${sample}-tmap.bam
  """
}
