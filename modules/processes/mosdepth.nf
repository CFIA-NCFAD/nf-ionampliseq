process MOSDEPTH_GENOME {
  tag "$sample"
  label 'process_low'

  conda 'bioconda::mosdepth=0.3.8'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0'
  } else {
    container 'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'
  }

  input:
  tuple val(sample), path(bam_bai), path(ref_fasta)

  output:
  tuple val(sample), path("*.per-base.bed.gz"), emit: bedgz
  path("*.global.dist.txt"), emit: mqc
  path("*.{txt,gz,csi,tsv}"), emit: results
  path("versions.yml"), emit: versions

  script:
  """
  mosdepth \\
      --fast-mode \\
      $sample \\
      ${bam_bai[0]}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
  END_VERSIONS
  """
}
