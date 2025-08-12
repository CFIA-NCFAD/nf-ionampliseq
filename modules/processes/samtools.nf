process SAMTOOLS_STATS {
  tag "$sample"
  
  conda 'bioconda::samtools=1.22'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/samtools:1.22--h96c455f_0'
  } else {
    container 'quay.io/biocontainers/samtools:1.22--h96c455f_0'
  }

  input:
  tuple val(sample), path(bam), path(ref_fasta)

  output:
  path('*.{flagstat,idxstats,stats}'), emit: stats
  path("versions.yml"), emit: versions

  script:
  """
  samtools flagstat *.bam > ${sample}.flagstat
  samtools idxstats *.bam > ${sample}.idxstats
  samtools stats *.bam > ${sample}.stats

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}

process SAMTOOLS_DEPTH {
  tag "$sample"
  
  conda 'bioconda::samtools=1.22'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/samtools:1.22--h96c455f_0'
  } else {
    container 'quay.io/biocontainers/samtools:1.22--h96c455f_0'
  }

  input:
  tuple val(sample), path(bam), path(ref_fasta)

  output:
  tuple val(sample), path('*-depths.tsv'), path(ref_fasta), emit: depth

  script:
  """
  samtools depth -aa -d 0 ${bam[0]} > ${sample}-depths.tsv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}

process SAMTOOLS_FASTQ {
  tag "$sample"

  label 'process_low'

  conda 'bioconda::samtools=1.22 conda-forge::pigz'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'
  }

  input:
  tuple val(sample), path(bam)

  output:
  tuple val(sample), path("*.fastq.gz"), emit: fastq_gz
  path("versions.yml"), emit: versions

  script:
  """
  samtools fastq $bam | pigz -c - > ${sample}.fastq.gz

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      pigz: \$(pigz --version 2>&1 | sed 's/pigz //g' )
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}