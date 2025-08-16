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
  echo "Major AF threshold: $major_allele_fraction"
  echo "Minor AF threshold: $minor_allele_fraction"
  echo "Using: AF, FRO, FAO, FDP (ignoring AO, RO, DP)"

  echo "Step 1: Analyzing multiallelic sites with Flow Evaluator data..."

  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AF]\\t[%FAO]\\t[%FRO]\\t[%FDP]\\n' "$vcf" \\
  | awk -v maj="$major_allele_fraction" -v min="$minor_allele_fraction" -v OFS='\\t' '
  {
    chrom=\$1; pos=\$2; ref=\$3; alt_str=\$4; af_str=\$5; fao_str=\$6; fro=\$7; fdp=\$8;
    
    # Parse alternates and Flow Evaluator metrics
    n_alts = split(alt_str, alts, ",");
    n_afs = split(af_str, afs, ",");
    split(fao_str, faos, ",");
    
    # Calculate total AF from Flow Evaluator
    total_af = 0;
    for(i = 1; i <= n_afs; i++) {
      total_af += (afs[i] + 0);
    }
    
    # Reference AF from Flow Evaluator
    if (fdp > 0) {
        ref_af = (fro + 0) / (fdp + 0);
    } else {
        ref_af = 0;
    }
    
    printf "# Site %s:%s - Total AF: %.3f, Ref AF: %.3f, FDP: %d\\n", chrom, pos, total_af, ref_af, fdp > "/dev/stderr";
    
    keep_alts = "";
    keep_afs = "";
    keep_faos = "";
    keep_indices = "";
    
    for(i = 1; i <= n_alts; i++) {
      current_af = (afs[i] + 0);
      current_alt = alts[i];
      current_fao = (faos[i] + 0);
      
      # Skip if below minimum threshold
      if(current_af < min) {
        printf "#   Alt %d (%s): AF=%.3f FAO=%d - BELOW_MIN_THRESHOLD\\n", i, current_alt, current_af, current_fao > "/dev/stderr";
        continue;
      }
      
      # Check if its an indel
      ref_len = length(ref);
      alt_len = length(current_alt);
      is_indel = (ref_len != alt_len) ? 1 : 0;
      
      # Include if: (AF >= min) AND (not indel OR AF >= major)
      # This excludes low-AF indels but keeps low-AF SNPs for ambiguity
      if(!is_indel || current_af >= maj) {
        if(keep_alts == "") {
          keep_alts = current_alt;
          keep_afs = current_af + "";
          keep_faos = current_fao + "";
          keep_indices = i + "";
        } else {
          keep_alts = keep_alts "," current_alt;
          keep_afs = keep_afs "," current_af;
          keep_faos = keep_faos "," current_fao;
          keep_indices = keep_indices "," i;
        }
        printf "#   Alt %d (%s): AF=%.3f FAO=%d - KEPT\\n", i, current_alt, current_af, current_fao > "/dev/stderr";
      } else {
        printf "#   Alt %d (%s): AF=%.3f FAO=%d - EXCLUDED (low AF indel)\\n", i, current_alt, current_af, current_fao > "/dev/stderr";
      }
    }
    
    # Determine final genotype based on kept alternates and total AF
    n_kept = split(keep_indices, kept_idx, ",");
    
    if(n_kept == 0) {
      # No alternates kept - reference
      final_gt = "0/0";
      action = "REF_ONLY";
      out_alts = ref;
      out_afs = ref_af + "";
    } else if(n_kept == 1) {
      # Single alternate
      single_af = (split(keep_afs, af_arr, ","), af_arr[1] + 0);
      if(total_af >= maj) {
          final_gt = "1/1";  # Homozygous alt (high total AF)
          action = "SINGLE_HOM";
      } else {
          final_gt = "0/1";  # Heterozygous with reference
          action = "SINGLE_HET";
      }
      out_alts = keep_alts;
      out_afs = keep_afs;
    } else if(n_kept == 2) {
      # Two alternates
      split(keep_afs, af_arr, ",");
      af1 = af_arr[1] + 0;
      af2 = af_arr[2] + 0;
      action = "DUAL_ALT_AMBIGUOUS";
      if(total_af >= maj) {
        final_gt = "1/2";  # Both alternates, creates ambiguity
      } else {
        # Both alternates meet minor threshold - keep both for ambiguity
        final_gt = "0/1";  # Heterozygous with reference, but ambiguous between two alternates
      }
      out_alts = keep_alts;
      out_afs = keep_afs;
    } else {
      # More than 2 alternates - keep alleles with AF > $minor_allele_fraction
      split(keep_alts, alt_arr, ",");
      split(keep_afs, af_arr, ",");
      
      # Sort by AF descending (simple bubble sort)
      for(i = 1; i <= n_kept; i++) {
        for(j = i+1; j <= n_kept; j++) {
          if((af_arr[j] + 0) > (af_arr[i] + 0)) {
            temp = af_arr[i]; af_arr[i] = af_arr[j]; af_arr[j] = temp;
            temp = alt_arr[i]; alt_arr[i] = alt_arr[j]; alt_arr[j] = temp;
          }
        }
      }
      
      out_alts = alt_arr[1] "," alt_arr[2];
      out_afs = af_arr[1] "," af_arr[2];
      final_gt = "1/2";  # Creates ambiguity code
      action = "MULTI_TO_DUAL";
    }
    
    printf "#   Decision: %s -> %s (action: %s)\\n", final_gt, out_alts, action > "/dev/stderr";
    
    # Output for further processing
    if(final_gt != "0/0" || total_af < min) {
      print chrom, pos, ref, out_alts, out_afs, final_gt, fdp, action;
    }
  }' > "${sample}.flow_eval_decisions.tsv"

  echo "Processed \$(wc -l < "${sample}.flow_eval_decisions.tsv") sites with Flow Evaluator metrics"

  echo "Step 2: Preparing annotations for bcftools..."

  # Create annotation file with proper columns for bcftools
  awk -v OFS='\\t' '{
      # CHROM, POS, REF, ALT, FLOW_AF, FLOW_GT, FLOW_ACTION, SELECTED
      print \$1, \$2, \$3, \$4, \$5, \$6, \$8, "1"
  }' "${sample}.flow_eval_decisions.tsv" > "${sample}.annotations.tsv"

  bgzip -c "${sample}.annotations.tsv" > "${sample}.annotations.tsv.gz"
  tabix -s1 -b2 -e2 "${sample}.annotations.tsv.gz"

  echo "Step 3: Applying Flow Evaluator based filtering..."

  # First split multiallelics to match our decisions
  bcftools norm -m -any "${vcf}" -Ou \\
  | bcftools annotate \\
      -a "${sample}.annotations.tsv.gz" \\
      -h <(printf '##INFO=<ID=FLOW_AF,Number=1,Type=Float,Description="Allele frequency from Flow Evaluator">\\n##INFO=<ID=FLOW_GT,Number=1,Type=String,Description="Genotype based on Flow Evaluator AF">\\n##INFO=<ID=FLOW_ACTION,Number=1,Type=String,Description="Action taken based on Flow Evaluator">\\n##INFO=<ID=SELECTED,Number=0,Type=Flag,Description="Selected for consensus">\\n') \\
      -c CHROM,POS,REF,ALT,INFO/FLOW_AF,INFO/FLOW_GT,INFO/FLOW_ACTION,INFO/SELECTED \\
      -Ou \\
  | bcftools view -i 'INFO/SELECTED=1' -Oz -o "${sample}.flow_selected.vcf.gz"

  bcftools index "${sample}.flow_selected.vcf.gz"

  echo "Step 4: Setting genotypes based on Flow Evaluator analysis..."

  # Use bcftools +setGT with correct syntax to set genotypes from FLOW_GT INFO field
  # Use temporary files to avoid corruption from overwriting compressed files
  
  # Set 1/1 genotypes where FLOW_GT=1/1
  bcftools +setGT "${sample}.flow_selected.vcf.gz" \\
      -Ov \\
      -o "${sample}.temp_1.vcf" \\
      -- -t q -n 'c:1/1' -i 'INFO/FLOW_GT="1/1"'

  # Set 0/1 genotypes where FLOW_GT=0/1  
  bcftools +setGT "${sample}.temp_1.vcf" \\
      -Ov \\
      -o "${sample}.temp_2.vcf" \\
      -- -t q -n 'c:0/1' -i 'INFO/FLOW_GT="0/1"'

  # Set 1/2 genotypes where FLOW_GT=1/2
  bcftools +setGT "${sample}.temp_2.vcf" \\
      -Ov \\
      -o "${sample}.temp_3.vcf" \\
      -- -t q -n 'c:1/2' -i 'INFO/FLOW_GT="1/2"'

  # Set 0/0 genotypes where FLOW_GT=0/0
  bcftools +setGT "${sample}.temp_3.vcf" \\
      -Ov \\
      -o "${sample}.flow_genotyped.vcf" \\
      -- -t q -n 'c:0/0' -i 'INFO/FLOW_GT="0/0"'

  echo "Step 5: Processing multiallelic sites for ambiguity codes..."

  # For sites with 1/2 genotype, we need to reconstruct the multiallelic record
  # This is complex with bcftools, so we'll create a special handling step

  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\tINFO/FLOW_ACTION\\n' "${sample}.flow_genotyped.vcf" \\
  | awk -v OFS='\\t' '\$5 == "1/2" {print \$1, \$2}' > "${sample}.multiallelic_sites.txt"

  if [ -s "${sample}.multiallelic_sites.txt" ]; then
      echo "Found \$(wc -l < "${sample}.multiallelic_sites.txt") true multiallelic sites for ambiguity handling"
      
      # For now, keep the 1/2 genotypes - bcftools consensus should handle these
      # by creating ambiguity codes (e.g., R for A/G, Y for C/T, etc.)
      cp "${sample}.flow_genotyped.vcf" "${sample}.final_processed.vcf"
  else
      echo "No true multiallelic sites found"
      cp "${sample}.flow_genotyped.vcf" "${sample}.final_processed.vcf"
  fi

  # Step 6: Optional frameshift filtering
  if [ "${filter_frameshift_variants}" = "true" ]; then
      echo "Step 6: Applying frameshift filtering..."
      
      bcftools filter "${sample}.final_processed.vcf" \\
          -e "TYPE != 'SNP' && (STRLEN(ALT) - STRLEN(REF)) % 3 != 0" \\
          -Ov -o "${sample}.bcftools_filt.vcf"
  else
      cp "${sample}.final_processed.vcf" "${sample}.bcftools_filt.vcf"
  fi

  # Step 7: Summary and examples
  echo "=== Processing Summary (Flow Evaluator Based) ==="
  echo "Input variants: \$(bcftools view -H "${vcf}" | wc -l)"
  echo "Final variants: \$(bcftools view -H "${sample}.bcftools_filt.vcf" | wc -l)"

  echo "Pre-Filter Genotype distribution:"
  bcftools query -f '[%GT]\\n' "${vcf}" | sort | uniq -c
  echo "Post-Filter Genotype distribution:"
  bcftools query -f '[%GT]\\n' "${sample}.bcftools_filt.vcf" | sort | uniq -c

  echo ""
  echo "Flow Evaluator Action distribution:"
  bcftools query -f '%INFO/FLOW_ACTION\\n' "${sample}.bcftools_filt.vcf" | sort | uniq -c

  echo "bcftools consensus -f reference.fasta -m low_coverage.bed ${sample}.bcftools_filt.vcf > ${sample}.consensus.fasta"

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

  output:
  tuple val(sample), path('*.bcftools.consensus.fasta'), emit: fasta
  path('*.merged.vcf.gz'), emit: vcf
  path("versions.yml"), emit: versions

  script:
  def consensus = "${sample}.bcftools.consensus.fasta"
  """
  # get low coverage depth mask BED file by filtering for regions with less than ${low_coverage}X
  zcat $mosdepth_per_base | awk '\$4<${low_coverage}' > low_cov.bed

  awk '/^>/ {print; next} {gsub(/[RYSWKMBDHVryswkmbdhv]/, "N"); print}' $fasta > ${fasta}.no_ambiguous.fasta

  # merge multiallelic sites into single record for bcftools consensus
  bcftools norm -m +any ${vcf} -Ou \\
  | bcftools +setGT -Ov -o ${vcf}.merged.vcf \\
    -- -t q -n 'c:0/1' -i 'INFO/FLOW_GT="0/1"'

  # Also force other genotypes to be preserved
  bcftools +setGT ${vcf}.merged.vcf -Ov -o ${vcf}.merged.vcf.tmp \\
    -- -t q -n 'c:1/1' -i 'INFO/FLOW_GT="1/1"'

  bcftools +setGT ${vcf}.merged.vcf.tmp -Ov -o ${vcf}.merged.vcf \\
    -- -t q -n 'c:1/2' -i 'INFO/FLOW_GT="1/2"'

  # Compress the final merged file
  bgzip -c ${vcf}.merged.vcf > ${vcf}.merged.vcf.gz
  tabix ${vcf}.merged.vcf.gz

  bcftools consensus \\
    -f ${fasta}.no_ambiguous.fasta \\
    -m low_cov.bed \\
    ${vcf}.merged.vcf.gz > $consensus

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
