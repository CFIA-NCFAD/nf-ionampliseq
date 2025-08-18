process QC_CSFV_FASTA {
  conda 'bioconda::biopython conda-forge::pandas conda-forge::rich conda-forge::typer conda-forge::edlib'
  // multicontainer for 
  // python=3.10,biopython=1.80,pandas=1.5.3,rich=12.6.0,typer=0.7.0,numpy=1.24.2,edlib=1.3.9
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  }

  input:
  path(fasta, stageAs: 'fastas/*')

  output:
  path("csfv_quality_summary.csv"), emit: qc_results
  path("csfv_quality_summary_mqc.txt"), emit: qc_results_mqc
  path('versions.yml'), emit: versions

  script:
  """
  qc_csfv_fasta.py \\
    fastas/ \\
    -v \\
    -o csfv_quality_summary.csv \\
    --output-mqc csfv_quality_summary_mqc.txt
  
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      qc_csfv_fasta.py: \$(qc_csfv_fasta.py --version)
  END_VERSIONS
  """
}
