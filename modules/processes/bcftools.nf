process BCFTOOLS_FILTER {
  tag "$sample"
  label 'process_low'

  conda 'bioconda::bcftools=1.20 conda-forge::gsl=2.7'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  }

  input:
  tuple val(sample), path(fasta), path(vcf)
  val(major_allele_fraction)
  val(minor_allele_fraction)
  val(filter_frameshift_variants)

  output:
  tuple val(sample), path(fasta), path(bcftools_filt_vcf), emit: vcf
  path("versions.yml"), emit: versions

  script:
  bcftools_filt_vcf = "${sample}.bcftools_filt.vcf"
  """
  bcftools norm \\
    --check-ref w \\
    -Ov \\
    -m- \\
    -f $fasta \\
    $vcf \\
    > norm.vcf

  # setGT plugin is used to set the genotype field (GT) based on the TVC Flow Evaluator calculated allele frequency (AF)
  bcftools +setGT \\
    norm.vcf \\
    -Ov \\
    -o setGT.major.vcf \\
    -- -t q -n 'c:1/1' -i 'FMT/AF >= ${major_allele_fraction}'

  bcftools +setGT \\
    setGT.major.vcf \\
    -Ov \\
    -o setGT.minor.vcf \\
    -- -t q  -n 'c:0/1' -i 'FMT/AF >= ${minor_allele_fraction} && FMT/AF < ${major_allele_fraction}'

  bcftools +setGT \\
    setGT.minor.vcf \\
    -Ov \\
    -o setGT.final.vcf \\
    -- -t q -n 'c:0/0' -i 'FMT/AF < ${minor_allele_fraction}'

  # Filter out frameshift variants if specified; may be too strict for Ion Torrent data since TVC should account for frameshifts in its error modeling.
  if [${filter_frameshift_variants} -eq true]; then
    bcftools filter \\
      setGT.final.vcf \\
      -e "TYPE != 'SNP' && ( (STRLEN(ALT) - STRLEN(REF)) % 3 ) != 0 
      || TYPE != 'SNP' && ( (STRLEN(ALT) - STRLEN(REF)) % 3 == 0
      && FMT/AF < ${minor_allele_fraction})" \\
      -Ov \\
      -o $bcftools_filt_vcf 
  else
    cp setGT.final.vcf $bcftools_filt_vcf 
  fi

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
  END_VERSIONS
  """
}

process BCFTOOLS_CONSENSUS {
  tag "$sample"
  label 'process_low'

  conda 'bioconda::bcftools=1.20 conda-forge::gsl=2.7'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  }

  input:
  tuple val(sample), path(fasta), path(vcf), path(mosdepth_per_base)
  val(low_coverage)
  val(major_allele_fraction)

  output:
  tuple val(sample), path('*.bcftools.consensus.fasta'), emit: fasta
  path("versions.yml"), emit: versions

  script:
  def consensus = "${sample}.bcftools.consensus.fasta"
  """
  # get low coverage depth mask BED file by filtering for regions with less than ${low_coverage}X
  zcat $mosdepth_per_base | awk '\$4<${low_coverage}' > low_cov.bed

  awk '/^>/ {print; next} {gsub(/[RYSWKMBDHVryswkmbdhv]/, "N"); print}' $fasta > ${fasta}.no_ambiguous.fasta

  bcftools filter \\
    -Oz \\
    -o no_low_af_indels.vcf.gz \\
    -e "TYPE != 'SNP' && FMT/AF < ${major_allele_fraction}" \\
    $vcf

  tabix no_low_af_indels.vcf.gz

  bcftools consensus \\
    -f ${fasta}.no_ambiguous.fasta \\
    -m low_cov.bed \\
    no_low_af_indels.vcf.gz > $consensus

  sed -i -E "s/^>(.*)/>${sample}/g" $consensus

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
  END_VERSIONS
  """
}

process BCFTOOLS_STATS {
  tag "$sample"
  label 'process_low'

  conda 'bioconda::bcftools=1.20 conda-forge::gsl=2.7'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  }
  input:
  tuple val(sample), path(fasta), path(vcf)

  output:
  path("*.bcftools_stats.txt"), emit: stats
  path("versions.yml"), emit: versions

  script:
  """
  bcftools stats -F $fasta $vcf > ${sample}.bcftools_stats.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
  END_VERSIONS
  """
}
