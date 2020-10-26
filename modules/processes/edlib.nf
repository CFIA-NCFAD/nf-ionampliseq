process EDLIB_ALIGN {
  tag "$sample"

  publishDir "${params.outdir}/consensus/edlib", 
             pattern: "*.txt", 
             mode: params.publish_dir_mode

  input:
  tuple val(sample), path(consensus_fasta), path(ref_fasta)

  output:
  tuple val(sample), path('*.txt'), emit: txt
  path '*.json', emit: json

  script:
  """
  edlib_align.py \\
    --sample $sample \\
    --consensus $consensus_fasta \\
    --reference $ref_fasta \\
    --output-aln ${sample}.edlib.txt \\
    --output-json ${sample}.json
  """
}

process EDLIB_MULTIQC {
  input:
  path(edlib_json)

  output:
  path('edlib_summary_mqc.txt')

  script:
  """
  edlib_multiqc.py ./ edlib_summary_mqc.txt
  """
}
