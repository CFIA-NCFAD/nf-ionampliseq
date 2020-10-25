process CONSENSUS {
  tag "$sample"
  publishDir "${params.outdir}/consensus", 
    pattern: "*.consensus.fasta",
    mode: params.publish_dir_mode

  input:
    tuple val(sample),
          path(vcf),
          path(ref_fasta),
          path(depths)
  output:
    tuple val(sample),
          path(vcf),
          path(ref_fasta),
          path(depths),
          path(consensus)

  script:
  consensus = "${sample}.consensus.fasta"
  """
  vcf_consensus_builder \\
    -v $vcf \\
    -d $depths \\
    -r $ref_fasta \\
    -o $consensus \\
    --low-coverage $params.low_coverage \\
    --no-coverage $params.no_coverage \\
    --low-cov-char $params.low_cov_char \\
    --no-cov-char $params.no_cov_char \\
    --sample-name $sample
  """
}
