process MASH_SCREEN {
  tag "$sample"
  publishDir "${params.outdir}/mash_screen",
             pattern: "*-mash_screen.tsv",
             mode: 'copy'
  publishDir "${params.outdir}/mash_screen",
             pattern: "*-top_ref.fasta",
             mode: 'copy'
  cpus 16
  
  input:
  tuple val(sample), path(reads), path(ref_fasta)

  output:
  tuple val(sample), path('*-mash_screen.tsv'), emit: results
  tuple val(sample), path('*-mash_screen-top_ref.fasta'), emit: top_ref

  script:
  """
  mash sketch -i -k ${params.mash_k} -s ${params.mash_s} $ref_fasta
  mash screen -p ${task.cpus} -w ${ref_fasta}.msh $reads | sort -grk1 > ${sample}-mash_screen.tsv
  cat ${sample}-mash_screen.tsv | head -n1 | awk '{print \$5}' > ${sample}-mash_screen-top_ref.txt
  seqtk subseq $ref_fasta ${sample}-mash_screen-top_ref.txt > ${sample}-mash_screen-top_ref.fasta
  """
}

process MASH_SCREEN_MULTIQC_SUMMARY {
  input:
  path(mash_screen_results)

  output:
  path('mash_screen_mqc.txt')

  script:
  """
  mash_screen_multiqc_summary.py ./ mash_screen_mqc.txt
  """
}