process SOFTWARE_VERSIONS {
  tag "Parse software version numbers"
  publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
      saveAs: { filename ->
                  if (filename.indexOf(".csv") > 0) filename
                  else null
              }

  output:
  path 'software_versions_mqc.yaml', emit: software_versions_yaml
  path "software_versions.csv"

  script:
  """
  echo $workflow.manifest.version > v_pipeline.txt
  echo $workflow.nextflow.version > v_nextflow.txt
  fastqc --version > v_fastqc.txt
  multiqc --version > v_multiqc.txt
  tmap --version 2> v_tmap.txt
  tvc --version > v_tvc.txt
  samtools --version > v_samtools.txt
  bcftools --version > v_bcftools.txt
  mash --version > v_mash.txt
  pigz --version 2> v_pigz.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}

process OUTPUT_DOCS {
  publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

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
