process SEQTK_SUBSEQ {
  tag "$sample"

  conda "bioconda::seqtk=1.5"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/seqtk:1.5--h577a1d6_1'
  } else {
    container 'quay.io/biocontainers/seqtk:1.5--h577a1d6_1'
  }

  input:
  tuple val(sample), path(ref_fasta), path(top_ref_txt)

  output:
  tuple val(sample), path('*.mash_screen-top_ref.fasta'), emit: top_ref_fasta
  path("versions.yml"), emit: versions

  script:
  """
  seqtk subseq $ref_fasta $top_ref_txt > ${sample}.mash_screen-top_ref.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      seqtk: \$(seqtk 2>&1 | grep -m1 '^Version:' | awk '{print \$2}')
  END_VERSIONS
  """
}
