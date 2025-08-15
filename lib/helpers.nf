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
  c_red = params.monochrome_logs ? '' : "\033[0;31m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  
  // Read schema JSON for parameter information
  def schemaFile = new File("${workflow.projectDir}/nextflow_schema.json")
  def schema = new groovy.json.JsonSlurper().parseText(schemaFile.text)
  
  log.info"""
  =${c_dim}=================================================================${c_reset}
  ${c_blue+c_bold}${workflow.manifest.name}${c_reset}   ~  version ${c_purple}${workflow.manifest.version}${c_reset}
  ${c_dim}==================================================================${c_reset}

    ${c_bold}Git info:${c_reset} $workflow.repository - $workflow.revision [$workflow.commitId]

  ${c_bold}Usage:${c_reset}

  The typical command for running the pipeline is as follows:
  
  \$ nextflow run ${workflow.manifest.name} \\
      ${c_red}--input '/path/to/iontorrent/*.bam'${c_reset} \\
      ${c_green}--outdir ${params.outdir}${c_reset} \\
      -profile docker # Recommended to run workflow with either Docker or Singularity enabled

  ${c_bold}Parameters:${c_reset}
"""

  // Generate help text from schema definitions
  schema.definitions.each { defName, defObj ->
    if (defName != "institutional_config_options" && defName != "max_job_request_options") {
      log.info "  ${c_bold}${defObj.title}:${c_reset}"
      log.info "    ${defObj.description}"
      log.info ""
      
      defObj.properties?.each { paramName, paramObj ->
        if (!paramObj.hidden) {
          def defaultValue = paramObj.default != null ? " (default: ${paramObj.default})" : ""
          log.info "    ${c_red}--${paramName}${c_reset}${defaultValue}"
          log.info "      ${paramObj.description}"
          if (paramObj.help_text) {
            // Convert markdown to terminal formatting
            def formattedHelpText = paramObj.help_text
              .replaceAll(/```bash\n([^`]+)```/, "${c_cyan}\$1${c_reset}")  // Code blocks in cyan
              .replaceAll(/```\n([^`]+)```/, "${c_cyan}\$1${c_reset}")      // Code blocks in cyan
              .replaceAll(/`([^`]+)`/, "${c_yellow}\$1${c_reset}")          // Inline code in yellow
              .replaceAll(/\*\*([^*]+)\*\*/, "${c_bold}\$1${c_reset}")     // Bold text
              .replaceAll(/\*([^*]+)\*/, "${c_dim}\$1${c_reset}")           // Italic text
            log.info "      ${formattedHelpText}"
          }
          log.info ""
        }
      }
    }
  }

  log.info """  ${c_bold}Other Options:${c_reset}
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

def onComplete(workflow, summary, custom_runName) {
    // Set up the e-mail variables
    def subject = "[CFIA-NCFAD/nf-ionampliseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[CFIA-NCFAD/nf-ionampliseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[CFIA-NCFAD/nf-ionampliseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[CFIA-NCFAD/nf-ionampliseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("${workflow.projectDir}/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("${workflow.projectDir}/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "${workflow.projectDir}", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("${workflow.projectDir}/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[CFIA-NCFAD/nf-ionampliseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[CFIA-NCFAD/nf-ionampliseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[CFIA-NCFAD/nf-ionampliseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[CFIA-NCFAD/nf-ionampliseq]${c_red} Pipeline completed with errors${c_reset}-"
    }
}
