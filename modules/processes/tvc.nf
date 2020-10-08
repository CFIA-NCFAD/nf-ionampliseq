process TVC {
  tag "$sample"

  publishDir "${params.outdir}/tvc/$sample",
             mode: 'copy'

  cpus 16
  memory 16.GB

  input:
  tuple val(sample),
        path(bam),
        path(ref_fasta),
        path(bed_file)

  output:
  tuple val(sample),
        path('*-tvc-postprocessed.{bam,bam.bai}'),
        path('tvc-outdir/*-small_variants.vcf'),
        path(ref_fasta), emit: vcf
  path 'tvc-outdir/'

  script:
  tvc_error_motifs_dir = (params.tvc_error_motifs_dir) ? "--error-motifs-dir ${params.tvc_error_motifs_dir}" : ""
  """
  cp $ref_fasta ref.fasta
  samtools faidx ref.fasta
  tvc \\
    --output-dir tvc-outdir \\
    --reference ref.fasta \\
    --input-bam ${bam[0]} \\
    --target-file $bed_file \\
    --trim-ampliseq-primers on \\
    --read-limit ${params.tvc_read_limit} \\
    --downsample-to-coverage ${params.tvc_downsample_to_coverage} \\
    --num-threads ${task.cpus} \\
    --postprocessed-bam ${sample}-tvc-postprocessed.bam \\
    --allow-complex \\
    --min-mapping-qv ${params.tvc_min_mapping_qv} \\
    --read-snp-limit ${params.tvc_read_snp_limit} \\
    --disable-filters \\
    $tvc_error_motifs_dir \\
    --output-vcf ${sample}-small_variants.vcf
  samtools index ${sample}-tvc-postprocessed.bam
  """
}
