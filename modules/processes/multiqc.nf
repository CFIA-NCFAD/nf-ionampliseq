process MULTIQC {
  label 'process_long'

  conda "bioconda::multiqc=1.30"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/multiqc:1.30--pyhdfd78af_0'
  } else {
    container 'quay.io/biocontainers/multiqc:1.30--pyhdfd78af_0'
  }

  input:
  path(multiqc_config)
  path(mqc_custom_config)
  path('fastqc/*')
  path('samtools/*')
  path('mosdepth/*')
  path('mash_screen/*')
  path('bcftools/*')
  path('edlib/*')
  path('consensus/*')
  path('blast/*')
  path('qc_csfv_fasta/*')
  path('software_versions/*')
  path(workflow_summary)

  output:
  path "*multiqc_report.html", emit: multiqc_report
  path "*_data"
  path "multiqc_plots"

  script:
  custom_runName = params.name
  if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
  }
  rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
  rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
  custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
  """
  multiqc --verbose -f $rtitle $rfilename $custom_config_file .
  """
}

process CONSENSUS_MULTIQC {
  conda 'bioconda::edlib bioconda::biopython conda-forge::pandas conda-forge::rich conda-forge::typer'
  // multicontainer for
  // python=3.10,biopython=1.80,pandas=1.5.3,rich=12.6.0,typer=0.7.0,numpy=1.24.2,edlib=1.3.9
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  }

  input:
  path(fastas)

  output:
  path('consensus_mqc.html')

  script:
  """
  consensus_fasta_multiqc.py $fastas consensus_mqc.html
  """
}
