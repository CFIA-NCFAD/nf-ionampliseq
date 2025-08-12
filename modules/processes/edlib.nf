process EDLIB_ALIGN {
  tag "$sample"
 
  conda 'bioconda::edlib bioconda::biopython conda-forge::pandas conda-forge::rich conda-forge::typer'
  // multicontainer for 
  // python=3.10,biopython=1.80,pandas=1.5.3,rich=12.6.0,typer=0.7.0,numpy=1.24.2,edlib=1.3.9
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  }

  input:
  tuple val(sample), path(consensus_fasta), path(ref_fasta)

  output:
  tuple val(sample), path('*.txt'), emit: txt
  path('*.json'), emit: json
  path('versions.yml'), emit: versions

  script:
  """
  edlib_align.py \\
    --sample $sample \\
    --consensus $consensus_fasta \\
    --reference $ref_fasta \\
    --output-aln ${sample}.edlib.txt \\
    --output-json ${sample}.json

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      edlib_align.py: \$(edlib_align.py --version)
  END_VERSIONS
  """
}

process EDLIB_MULTIQC {
  conda 'bioconda::edlib bioconda::biopython conda-forge::pandas conda-forge::rich conda-forge::typer'
  // multicontainer for 
  // python=3.10,biopython=1.80,pandas=1.5.3,rich=12.6.0,typer=0.7.0,numpy=1.24.2,edlib=1.3.9
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-4627b14ef9fdc900fafc860cd08c19d7fbb92f43:3aae397a063ad1c6f510775118a2c5e6acbc027d-0'
  }
  
  input:
  path(edlib_json)

  output:
  path('edlib_summary_mqc.txt')

  script:
  """
  edlib_multiqc.py ./ edlib_summary_mqc.txt
  """
}
