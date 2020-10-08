process FASTQC {
  tag "$sample"
  label 'process_medium'
  publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
      saveAs: { filename ->
                    filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
              }

  input:
  tuple val(sample), path(reads)

  output:
  path "*_fastqc.{zip,html}"

  script:
  """
  fastqc --quiet --threads $task.cpus $reads
  """
}
