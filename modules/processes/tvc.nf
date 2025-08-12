process TVC {
  tag "$sample"
  label 'process_medium'

  container 'ghcr.io/cfia-ncfad/nf-ionampliseq:2.0.0'

  input:
  tuple val(sample), path(bam), path(ref_fasta, stageAs: 'fasta/*'), path(bed_file)
  path(tvc_error_motifs_dir)
  val(trim_primers)

  output:
  tuple val(sample), path('*-tvc-postprocessed.{bam,bam.bai}'), path('tvc-outdir/*-small_variants.vcf'), path(ref_fasta), emit: vcf
  path 'tvc-outdir/', emit: outdir
  path("versions.yml"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def primer_trim_args = trim_primers ? "--trim-ampliseq-primers on --target-file $bed_file" : ""
  """
  cp $ref_fasta ref.fasta
  samtools faidx ref.fasta
  tvc \\
    $args \\
    $primer_trim_args \\
    --output-dir tvc-outdir \\
    --reference ref.fasta \\
    --input-bam ${bam[0]} \\
    --read-limit ${params.tvc_read_limit} \\
    --downsample-to-coverage ${params.tvc_downsample_to_coverage} \\
    --num-threads ${task.cpus} \\
    --postprocessed-bam ${sample}-tvc-postprocessed.bam \\
    --min-mapping-qv ${params.tvc_min_mapping_qv} \\
    --read-snp-limit ${params.tvc_read_snp_limit} \\
    --error-motifs-dir $tvc_error_motifs_dir \\
    --output-vcf ${sample}-small_variants.vcf
  samtools index ${sample}-tvc-postprocessed.bam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    tvc: \$(tvc --version 2>&1 | head -n1 | sed 's/^tvc //; s/ .*//')
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
