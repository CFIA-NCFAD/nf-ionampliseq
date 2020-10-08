process MOSDEPTH_GENOME {
  tag "$sample"
  label 'process_medium'
  publishDir "${params.outdir}/mosdepth", 
             mode: 'copy',
             saveAs: { filename ->
               if (filename.endsWith(".pdf")) "plots/$filename"
               else if (filename.endsWith(".tsv")) "plots/$filename"
               else filename
             }

  input:
  tuple val(sample), path(bam), path(ref_fasta)

  output:
  path "*.global.dist.txt", emit: mqc
  path "*.{txt,gz,csi,tsv}"

  script:
  """
  mosdepth \\
      --by 200 \\
      --fast-mode \\
      $sample \\
      ${bam[0]}
  """
}
