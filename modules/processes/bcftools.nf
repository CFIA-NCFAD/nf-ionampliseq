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
  def bcftools_frameshift_filter = filter_frameshift_variants ? "TYPE != 'SNP' && ( (STRLEN(ALT) - STRLEN(REF)) % 3 ) != 0 || TYPE != 'SNP' && ( (STRLEN(ALT) - STRLEN(REF)) % 3 == 0 && " : ""
  """
  bcftools norm \\
    --check-ref w \\
    -Ov \\
    -m- \\
    -f $fasta \\
    $vcf \\
    > norm.vcf

  # TSV with CHROM,POS,REF,ALT,AD
  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%RO,%AO]\\n' norm.vcf > ad.tsv

  # header file defining the AD FORMAT field
  cat <<-END_AD_HEADER > ad.hdr 
  ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
  END_AD_HEADER

  # bgzip and tabix index the TSV
  # Must be CHROM, POS, REF, ALT, AD (no FORMAT/ prefix here)
  bgzip -c ad.tsv > ad.tsv.gz
  tabix -s1 -b2 -e2 ad.tsv.gz

  # annotate the VCF
  bcftools annotate \
    -a ad.tsv.gz \
    -c CHROM,POS,REF,ALT,FMT/AD \
    -h ad.hdr \
    norm.vcf \
    > with_ad.vcf

  bcftools +fill-tags \\
    with_ad.vcf \\
    -Ov \\
    -o filled.vcf \\
    -- -t all

  bcftools +setGT \\
    filled.vcf \\
    -Ov \\
    -o setGT.major.vcf \\
    -- -t q -n 'c:1/1' -i 'FMT/VAF >= ${major_allele_fraction}'

  bcftools +setGT \\
    setGT.major.vcf \\
    -Ov \\
    -o setGT.minor.vcf \\
    -- -t q  -n 'c:0/1' -i 'FMT/VAF >= ${minor_allele_fraction} && FMT/VAF < ${major_allele_fraction}'

  bcftools +setGT \\
    setGT.minor.vcf \\
    -Ov \\
    -o setGT.final.vcf \\
    -- -t q -n 'c:0/0' -i 'FMT/VAF < ${minor_allele_fraction}'

  bcftools filter \\
    setGT.final.vcf \\
    -e "${bcftools_frameshift_filter} FMT/VAF < ${minor_allele_fraction})" \\
    -Ov \\
    -o $bcftools_filt_vcf 

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
    -e "TYPE != 'SNP' && FMT/VAF < ${major_allele_fraction}" \\
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
