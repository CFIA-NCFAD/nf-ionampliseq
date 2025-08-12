process OUTPUT_DOCS {
  conda 'bioconda::multiqc=1.30'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/multiqc:1.30--pyhdfd78af_0'
  } else {
    container 'quay.io/biocontainers/multiqc:1.30--pyhdfd78af_0'
  }

  input:
  path output_docs
  path images

  output:
  path "results_description.html"

  script:
  """
  markdown_to_html.py $output_docs -o results_description.html
  """
}


//
// Generate workflow version string
//
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
}

//
// Get software versions for pipeline
//
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file)
    return yaml.dumpAsMap(versions).trim()
}

//
// Get workflow version for pipeline
//
def workflowVersionToYAML() {
    return """
    Workflow:
        ${workflow.manifest.name}: ${getWorkflowVersion()}
        Nextflow: ${workflow.nextflow.version}
    """.stripIndent().trim()
}

//
// Get channel of software versions used in pipeline in YAML format
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions.unique().map { version -> 
      if (version == null) {
        return null
      } else {
        return processVersionsFromYAML(version)
      }
    }
    .filter { it != null }
    .unique().mix(Channel.of(workflowVersionToYAML()))
}

//
// Get software versions for MultiQC
//
def softwareVersionsMultiqc(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions_map = yaml.load(yaml_file)
    versions_map["Workflow"] = [
      "${workflow.manifest.name}": "${getWorkflowVersion()}",
      "Nextflow": "${workflow.nextflow.version}"
    ]
    def versions_section = ''
    
    // Start the HTML table with CSS styling
    versions_section += "    <style>\n"
    versions_section += "    #nf-core-versions tr:nth-child(even) {\n"
    versions_section += "        background-color: #f2f2f2;\n"
    versions_section += "    }\n"
    versions_section += "    </style>\n"
    versions_section += "    <table class=\"table\" style=\"width:100%\" id=\"nf-core-versions\">\n"
    versions_section += "        <thead>\n"
    versions_section += "            <tr>\n"
    versions_section += "                <th> Process Name </th>\n"
    versions_section += "                <th> Software </th>\n"
    versions_section += "                <th> Version  </th>\n"
    versions_section += "            </tr>\n"
    versions_section += "        </thead>\n"
    versions_section += "    <tbody>\n"
    
    // Process each entry in the versions map
    versions_map
        .keySet()
        .sort()
        .each { process_name ->
            def process_versions = versions_map.get(process_name)
            if (process_versions) {
                process_versions
                    .keySet()
                    .sort()
                    .each { software ->
                        def version = process_versions.get(software) ?: 'N/A'
                        versions_section += "    <tr>\n"
                        versions_section += "        <td><samp>${process_name}</samp></td>\n"
                        versions_section += "        <td><samp>${software}</samp></td>\n"
                        versions_section += "        <td><samp>${version}</samp></td>\n"
                        versions_section += "    </tr>\n"
                    }
            }
        }
    
    // Close the table
    versions_section += "    </tbody>\n"
    versions_section += "    </table>\n"
    
    // Build the YAML content
    def yaml_file_text = "plot_type: 'html'\n" as String
    yaml_file_text     += "data: |-\n"
    yaml_file_text     += "${versions_section}"
    yaml_file_text     += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    yaml_file_text     += "id: 'software_versions-module'\n"
    yaml_file_text     += "section_name: '${workflow.manifest.name} Software Versions'\n"
    yaml_file_text     += "description: 'Software versions are collected at run time from the software output.'\n"
    
    return yaml_file_text
}

//
// Get workflow summary for MultiQC
//
def paramsSummaryMultiqc(summary_params) {
    def summary_section = ''
    summary_params
        .keySet()
        .each { group ->
            def group_params = summary_params.get(group)
            // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>${group}</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                group_params
                    .keySet()
                    .sort()
                    .each { param ->
                        summary_section += "        <dt>${param}</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                    }
                summary_section += "    </dl>\n"
            }
        }

    def yaml_file_text = "id: '${workflow.manifest.name.replace('/', '-')}-summary'\n" as String
    yaml_file_text     += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_file_text     += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    yaml_file_text     += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    yaml_file_text     += "plot_type: 'html'\n"
    yaml_file_text     += "data: |\n"
    yaml_file_text     += "${summary_section}"

    return yaml_file_text
}
