process MASH_SKETCH {
  label 'process_low'
  conda 'bioconda::mash=2.3'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mash:2.3--hb105d93_9'
  } else {
    container 'quay.io/biocontainers/mash:2.3--hb105d93_9'
  }

  input:
  path(ref_fasta)

  output:
  tuple path(ref_fasta), path('*.msh'), emit: sketch
  path("versions.yml"), emit: versions

  script:
  """
  mash sketch -i -k ${params.mash_k} -s ${params.mash_s} $ref_fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      mash: \$(mash --version 2>&1)
  END_VERSIONS
  """
}

process MASH_SCREEN {
  tag "$sample"
  label 'process_low'
  conda 'bioconda::mash=2.3'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mash:2.3--hb105d93_9'
  } else {
    container 'quay.io/biocontainers/mash:2.3--hb105d93_9'
  }

  input:
  tuple val(sample), path(reads), path(ref_fasta), path(mash_sketch)

  output:
  tuple val(sample), path('*-mash_screen.tsv'), emit: results
  tuple val(sample), path('*-mash_screen-top_ref.txt'), emit: top_ref_txt
  path("versions.yml"), emit: versions

  script:
  """

  mash screen -p ${task.cpus} -w ${ref_fasta}.msh $reads | sort -grk1 > ${sample}-mash_screen.tsv
  cat ${sample}-mash_screen.tsv | head -n1 | awk '{print \$5}' > ${sample}-mash_screen-top_ref.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      mash: \$(mash --version 2>&1)
  END_VERSIONS
  """
}

process MASH_SCREEN_MULTIQC_SUMMARY {
  conda 'bioconda::edlib bioconda::biopython conda-forge::pandas conda-forge::rich conda-forge::typer'
  // multicontainer for
  // python=3.10,biopython=1.80,pandas=1.5.3,rich=12.6.0,typer=0.7.0,numpy=1.24.2,edlib=1.3.9
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  }

  input:
  path(mash_screen_results)

  output:
  path('mash_screen_mqc.txt')

  script:
  """
  mash_screen_multiqc_summary.py ./ mash_screen_mqc.txt
  """
}
