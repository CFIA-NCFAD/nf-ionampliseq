#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { 
  helpMessage; 
  checkHostname; 
  check_rundir; 
  check_sample_sheet 
  } from './lib/helpers'

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// If using Slurm profile then check that Slurm queue is specified, otherwise, exit with error
if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  exit 1, "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
}

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['References multi FASTA']           = params.ref_fasta
summary['AmpliSeq BED file']                = params.bed_file
summary['Mash kmer length']                 = params.mash_k
summary['Mash # of hashes']                 = params.mash_s
summary['TVC error motifs directory']       = params.tvc_error_motifs_dir
summary['TVC downsample to X coverage']     = params.tvc_downsample_to_coverage
summary['TVC min mapping Q value']          = params.tvc_min_mapping_qv
summary['TVC read SNP limit']               = params.tvc_read_snp_limit
summary['TVC read limit']                   = params.tvc_read_limit
summary['Consensus low coverage']           = params.low_coverage
summary['Consensus low coverage character'] = params.low_cov_char
summary['Consensus no coverage']            = params.no_coverage
summary['Consensus no coverage character']  = params.no_cov_char
summary['Sample sheet']                     = params.sample_sheet
summary['Max Resources']                    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) {
    summary['Container'] = "$workflow.containerEngine - $workflow.container"
}
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) {
    summary['Config Profile Description'] = params.config_profile_description
}
if (params.config_profile_contact) {
    summary['Config Profile Contact']     = params.config_profile_contact
}
if (params.config_profile_url) {
    summary['Config Profile URL']         = params.config_profile_url
}
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

// Include processes
include { OUTPUT_DOCS; SOFTWARE_VERSIONS } from './modules/processes/docs'
include { FASTQC } from './modules/processes/fastqc'
include { MULTIQC } from './modules/processes/multiqc'
include { CHECK_SAMPLE_SHEET; BAM_TO_FASTQ; FILTER_BED_FILE } from './modules/processes/misc'
include { MASH_SCREEN; MASH_SCREEN_MULTIQC_SUMMARY } from './modules/processes/mash'
include { TMAP } from './modules/processes/tmap'
include { SAMTOOLS_DEPTH; SAMTOOLS_STATS } from './modules/processes/samtools'
include { TVC } from './modules/processes/tvc'
include { BCFTOOLS_VCF_FILTER; BCFTOOLS_VCF_NORM } from './modules/processes/bcftools'
include { MOSDEPTH_GENOME } from './modules/processes/mosdepth'
include { CONSENSUS } from './modules/processes/consensus'
include { COVERAGE_PLOT } from './modules/processes/plotting'


workflow {
  if (params.sample_sheet) {
    sample_sheet_file = file(params.sample_sheet, checkIfExists: true)
    ch_input = Channel.from(sample_sheet_file) \
      | CHECK_SAMPLE_SHEET \
      | splitCsv(header: ['sample', 'bam'], sep: ',', skip: 1) \
      | map { check_sample_sheet(it) }
    if (params.panel && params.panels && !params.panels.containsKey(params.panel)) {
      exit 1, "The specified AmpliSeq panel '$params.panel' is not a valid panel. You must specify one of the following: ${params.panels.keySet().join(', ')}"
    }
    if (params.panel == null && (params.ref_fasta == null || params.bed_file == null)) {
      exit 1, "You must either specify a built-in AmpliSeq panel ('--panel' one of: ${params.panels.keySet().join(', ')}) or specify a valid reference FASTA file ('--ref_fasta') and BED file ('--bed_file') for an AmpliSeq panel."
    }
    if (params.panel && params.panels && params.panels.containsKey(params.panel)) {
      panel = params.panels.get(params.panel)
      ch_ref_fasta = Channel.from(file(panel.fasta, checkIfExists: true))
      ch_bed_file = Channel.from(file(panel.bed_file, checkIfExists: true))
    }
    if (params.ref_fasta && params.bed_file) {
      ch_ref_fasta = Channel.from(file(params.ref_fasta, checkIfExists: true))
      ch_bed_file = Channel.from(file(params.bed_file, checkIfExists: true))
    }
  } else if (params.rundir) {
    rundir = file(params.rundir, checkIfExists: true)
    (samples, ref_panel) = check_rundir(rundir)
    if (ref_panel == 'CSFV_AmpliSeq') {
      panel = params.panels.get('csf')
      log.info "Using AmpliSeq panel 'csf' given that reference was found to be '$ref_panel' in rundir '$rundir'. Panel='$panel'"
    } else if (ref_panel ~== 'FMDV_AmpliSeq') {
      panel = params.panels.get('fmd')
    } else {
      exit 1, "Not sure what AmpliSeq panel to use given reference '$ref_panel' found in 'ion_params_00.json'"
    }
    ch_ref_fasta = Channel.from(file(panel.fasta, checkIfExists: true))
    ch_bed_file = Channel.from(file(panel.bed_file, checkIfExists: true))
    ch_input = Channel.from(samples)
  } else { 
    exit 1, "Sample sheet tab-delimited file not specified! Please specify path to sample_sheet with '--sample_sheet /path/to/sample_sheet.tsv'"
  }

  ch_input | BAM_TO_FASTQ \
           | combine(ch_ref_fasta) \
           | MASH_SCREEN
  ch_input | join(MASH_SCREEN.out.top_ref) \
           | TMAP \
           | (SAMTOOLS_STATS & MOSDEPTH_GENOME & SAMTOOLS_DEPTH)
  MASH_SCREEN.out.results | combine(ch_bed_file) | FILTER_BED_FILE
  MASH_SCREEN.out.results \
    | collect \
    | map { [ it[1] ] } \
    | MASH_SCREEN_MULTIQC_SUMMARY
  BAM_TO_FASTQ.out | FASTQC

  TMAP.out | join(FILTER_BED_FILE.out) | TVC
  TVC.out.vcf | BCFTOOLS_VCF_NORM | BCFTOOLS_VCF_FILTER
  BCFTOOLS_VCF_FILTER.out \
    | join(SAMTOOLS_DEPTH.out) \
    | map { [it[0], it[2], it[3], it[4]]} \
    | (CONSENSUS & COVERAGE_PLOT)

  // Stage config files
  ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
  ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
  ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
  ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

  Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt style=\"width:240px !important;\">$k</dt><dd style=\"margin-left:260px !important;\"><samp>${v != null ? v : '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'peterk87-nf-ionampliseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'peterk87/nf-ionampliseq Workflow Summary'
    section_href: 'https://github.com/peterk87/nf-ionampliseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

  SOFTWARE_VERSIONS()
  OUTPUT_DOCS(ch_output_docs, ch_output_docs_images)
  MULTIQC(
    ch_multiqc_config,
    ch_multiqc_custom_config.collect().ifEmpty([]),
    FASTQC.out.collect().ifEmpty([]),
    SAMTOOLS_STATS.out.collect().ifEmpty([]),
    MOSDEPTH_GENOME.out.mqc.collect().ifEmpty([]),
    MASH_SCREEN_MULTIQC_SUMMARY.out.collect(),
    SOFTWARE_VERSIONS.out.software_versions_yaml.collect(),
    ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
  )
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[peterk87/nf-ionampliseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[peterk87/nf-ionampliseq] FAILED: $workflow.runName"
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
                log.warn "[peterk87/nf-ionampliseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[peterk87/nf-ionampliseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[peterk87/nf-ionampliseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[peterk87/nf-ionampliseq] Sent summary e-mail to $email_address (mail)"
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
        log.info "-${c_purple}[peterk87/nf-ionampliseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[peterk87/nf-ionampliseq]${c_red} Pipeline completed with errors${c_reset}-"
    }

}
