/*
* Utility processes
*/

process SAMPLE_INFO_FROM_BAM {
  tag "$bam"

  conda 'bioconda::pysam bioconda::biopython conda-forge::pandas conda-forge::rich conda-forge::typer'
  // biocontainers/multi-package-containers image for 
  // python=3.9,pysam,biopython,click,pandas,numpy,ete3
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-c31f862d526d0a07196020b22083c45f8ddb4f0d:5d34a7705065087f5129c939eea53217c22e38df-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-c31f862d526d0a07196020b22083c45f8ddb4f0d:5d34a7705065087f5129c939eea53217c22e38df-0'
  }


  input:
  path(bam)

  output:
  tuple path(bam), path('sample_name.txt'), path('ampliseq_panel.txt'), emit: sample_info
  path('*.tsv'), emit: tsv
  path('versions.yml'), emit: versions

  script:
  """
  sample_info_from_bam.py \\
    -i $bam \\
    -o sample_name.txt \\
    -p ampliseq_panel.txt \\
    --write-sample-info

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      sample_info_from_bam.py: \$(sample_info_from_bam.py --version)
  END_VERSIONS
  """
}

process CHECK_SAMPLE_SHEET {
  tag "$samplesheet"
  
  conda 'bioconda::edlib bioconda::biopython conda-forge::pandas conda-forge::rich conda-forge::typer'
  // multicontainer for 
  // python=3.10,biopython=1.80,pandas=1.5.3,rich=12.6.0,typer=0.7.0,numpy=1.24.2,edlib=1.3.9
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  }

  input:
  path(samplesheet)

  output:
  path "samplesheet_reformat.csv", emit: samplesheet
  path('versions.yml'), emit: versions

  script:
  """
  check_sample_sheet.py $samplesheet samplesheet_reformat.csv
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      check_sample_sheet.py: \$(check_sample_sheet.py --version)
  END_VERSIONS
  """
}

process FILTER_BED_FILE {
  tag "$sample"

  input:
  tuple val(sample), path(top_ref), path(bed_file)

  output:
  tuple val(sample), path("*-filtered.bed")

  script:
  """
  REF=\$(head -n1 $top_ref | awk '{print \$5}')
  head -n1 $bed_file > ${sample}-filtered.bed
  grep "^\$REF\\b" $bed_file >> ${sample}-filtered.bed
  """
}
