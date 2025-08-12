// modified from nf-core/modules
// https://github.com/nf-core/modules/blob/ab0581ef6b11c5e62c641c6194a847fb622987fa/modules/nf-core/fastqc/main.nf#L1
process FASTQC {
  tag "$sample"
  
  conda 'bioconda::fastqc=0.12.1'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'
  } else {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
  }
  
  input:
  tuple val(sample), path(reads)

  output:
  path("*_fastqc.{zip,html}"), emit: results
  path('versions.yml'), emit: versions

  script:
  def args = task.ext.args ?: ''
  // The total amount of allocated RAM by FastQC is equal to the number of threads defined (--threads) time the amount of RAM defined (--memory)
  // https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/fastqc#L211-L222
  // Dividing the task.memory by task.cpu allows to stick to requested amount of RAM in the label
  def memory_in_mb = task.memory ? task.memory.toUnit('MB') / task.cpus : null
  // FastQC memory value allowed range (100 - 10000)
  def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)
  """
  fastqc \\
    ${args} \\
    --threads ${task.cpus} \\
    --memory ${fastqc_memory} \\
    $reads

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    fastqc: \$(fastqc --version | sed '/FastQC v/!d; s/.*v//')
  END_VERSIONS
  """
}
