process BLASTN {
  tag "$sample"
  label 'process_medium'
  
  conda 'bioconda::blast=2.16.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/blast:2.16.0--hc155240_3'
  } else {
    container 'quay.io/biocontainers/blast:2.16.0--hc155240_3'
  }

  input:
  tuple val(sample), path(fasta), val(db_prefix), path(db)

  output:
  tuple val(sample), path("*.blast.tsv"), emit: tsv
  path("versions.yml"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  """
  blastn \\
    ${args} \\
    -num_threads ${task.cpus} \\
    -query "${fasta}" \\
    -db "${db}/${db_prefix}" \\
    -out "${sample}.blast.tsv" \\

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    blastn: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
  END_VERSIONS
  """
}

process BLASTN_COVERAGE {
  label 'process_medium'
  conda 'bioconda::biopython conda-forge::pandas conda-forge::rich conda-forge::typer'
  // multicontainer for 
  // python=3.10,biopython=1.80,pandas=1.5.3,rich=12.6.0,typer=0.7.0,numpy=1.24.2,edlib=1.3.9
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  }

  input:
  path(blast_tsv, stageAs: 'blast_tsvs/*')

  output:
  path("blast_coverage_summary.tsv"), emit: tsv
  path("blast_coverage_mqc.txt"), emit: mqc_table
  path("versions.yml"), emit: versions

  script:
  """
  blast_coverage.py blast_tsvs/ -o .

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    blast_coverage: \$(blast_coverage.py --version 2>&1 | sed 's/^.*blast_coverage.py: //; s/ .*\$//')
  END_VERSIONS
  """
}
