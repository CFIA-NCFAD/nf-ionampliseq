/*
 * -------------------------------------------------
 *  peterk87/nf-ionampliseq helper functions
 * -------------------------------------------------
 */

def helpMessage() {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_bold = params.monochrome_logs ? '' : "\033[1m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_block = params.monochrome_logs ? '' : "\033[3m";
  c_ul = params.monochrome_logs ? '' : "\033[4m";
  c_black = params.monochrome_logs ? '' : "\033[0;30m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_bul = params.monochrome_logs ? '' : c_bold + c_ul;
  log.info"""
  =${c_dim}=================================================================${c_reset}
  ${c_blue+c_bold}${workflow.manifest.name}${c_reset}   ~  version ${c_purple}${workflow.manifest.version}${c_reset}
  ${c_dim}==================================================================${c_reset}

    ${c_ul}Git info:${c_reset} $workflow.repository - $workflow.revision [$workflow.commitId]

  ${c_bul}Usage:${c_reset}

  The typical command for running the pipeline is as follows:
  
  \$ nextflow run ${workflow.manifest.name} \\
      ${c_red}--sample_sheet "${params.sample_sheet}"${c_reset} \\
      ${c_green}--outdir ${params.outdir}${c_reset} \\
      -profile docker # Recommended to run workflow with either Docker or Singularity enabled

  ${c_bul}Mandatory Options:${c_reset}
    ${c_red}--sample_sheet${c_reset}   Sample sheet CSV file (default: ${c_red}"${params.sample_sheet}"${c_reset})

  Mash Screen Options:
    --mash_k          Mash sketch kmer size (default: ${params.mash_k})
    --mash_s          Mash sketch number of sketch hashes (default: ${params.mash_s})
  
  TVC Options:
    --tvc_error_motifs_dir        Directory with Ion Torrent TVC error motifs (default: ${params.tvc_error_motifs_dir})
    --tvc_read_limit              TVC read limit (default: $params.tvc_read_limit)
    --tvc_downsample_to_coverage  TVC downsample to at most X coverage (default: X=$params.tvc_downsample_to_coverage)
    --tvc_min_mapping_qv          TVC min mapping quality value (default: $params.tvc_min_mapping_qv)
    --tvc_read_snp_limit          TVC: do not use reads with number of SNPs about this (default: $params.tvc_read_snp_limit)

  ${c_bul}Cluster Options:${c_reset}
    --slurm_queue     Name of SLURM queue to run workflow on; use with ${c_dim}-profile slurm${c_reset}

  ${c_bul}Other Options:${c_reset}
    ${c_green}--outdir${c_reset}          The output directory where the results will be saved
                      (default: ${c_green}${params.outdir}${c_reset})
    -w/--work-dir     The temporary directory where intermediate data will be 
                      saved (default: ${workflow.workDir})
    -profile          Configuration profile to use. [standard, singularity, 
                      conda, slurm] (default '${workflow.profile}')
    --tracedir        Pipeline run info output directory (default: 
                      ${params.tracedir})
  """.stripIndent()
}

def checkHostname() {
  def c_reset = params.monochrome_logs ? '' : "\033[0m"
  def c_white = params.monochrome_logs ? '' : "\033[0;37m"
  def c_red = params.monochrome_logs ? '' : "\033[1;91m"
  def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
  if (params.hostnames) {
    def hostname = "hostname".execute().text.trim()
    params.hostnames.each { prof, hnames ->
      hnames.each { hname ->
        if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
          log.error "====================================================\n" +
                  "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                  "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                  "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                  "============================================================"
        }
      }
    }
  }
}

def check_sample_sheet(LinkedHashMap sample_sheet) {
  // Check that each entry from a sample sheet has a BAM file that exists
  bam = sample_sheet.bam ? file(sample_sheet.bam, checkIfExists: true) : null
  if (bam == null) {
    exit 1, "BAM file not specified for ${sample_sheet.sample} does not exist! Please check the sample sheet '$params.sample_sheet'"
  }
  return [ sample_sheet.sample, bam ]
}

def check_rundir(rundir) {
  // Check an Ion Torrent run directory for BAM files and an ion_params_00.json file
  // Parse the ion_params_00.json file for sample names and the appropriate 
  // AmpliSeq panel to use for analysis.
  ion_params_file = null
  bam_map = [:]
  rundir.eachFile {
    fname = it.getName()
    m = (fname =~ /(IonCode_\d+)_rawlib\.bam$/)
    if (m) {
      bam_map[m.group(1)] = file(it)
    }
    if (fname ==~ /ion_params_00\.json$/) {
      ion_params_file = file(it)
    }
  }
  if (ion_params_file == null) {
    exit 1, "Cannot find the 'ion_params_00.json' file in ${rundir}. Please verify that you have specified an exported Ion Torrent sequencing run base directory."
  } else {
    log.info "Found '$ion_params_file' in '$rundir'"
  }
  if (bam_map.size() == 0) {
    exit 1, "Could not find any BAM files matching filename pattern 'IonCode_*_rawlib.bam' in $rundir. Please verify that you have specified an exported Ion Torrent sequencing run base directory."
  } else {
    log.info "Found ${bam_map.size()} BAM files in '$rundir'"
  }
  slurper = new groovy.json.JsonSlurper()
  ion_params_text = ion_params_file.text
  result  = slurper.parseText(ion_params_file.text)
  log.info "Parsed JSON $ion_params_file (length=${ion_params_text.size()})"
  sample_to_barcode = result.experimentAnalysisSettings.barcodedSamples
  samples = []
  reference = [:]
  sample_to_barcode.each { sample, value ->
    bc = value.barcodes[0]
    sample_info = value.barcodeSampleInfo.get(bc)
    r = sample_info.reference
    if (reference.containsKey(r)) {
      reference[r] << sample
    } else {
      reference[r] = [sample]
    }
    
    bam = bam_map.get(bc)
    if (bam == null) {
      exit 1, "Could not find BAM file for sample '$sample' with barcode '$bc'. Please verify that your Ion Torrent sequencing run was exported correctly."
    }
    samples << [sample, bam]
  }
  ref_panel = reference.keySet()[0]
  if (reference.get(ref_panel).size() != samples.size()) {
    exit 1, "Not all ${samples.size()} samples are the same reference '$ref_panel'!"
  }
  log.info "Reference determined to be '${ref_panel}'"
  log.info "Found ${samples.size()} barcoded samples"
  return [samples, ref_panel]
}
