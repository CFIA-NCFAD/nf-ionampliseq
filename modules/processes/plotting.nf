process COVERAGE_PLOT {
  tag "$sample"
  publishDir "${params.outdir}/plots", 
    pattern: '*.pdf',
    mode: 'copy'

  input:
  tuple val(sample),
        path(vcf),
        path(ref_fasta),
        path(depths)
  
  output:
  path("*.pdf")

  script:
  plot_base_filename = "coverage_plot-${sample}"
  """
  plot_coverage.py --sample-name $sample -d $depths -o ${plot_base_filename}.pdf
  plot_coverage.py --sample-name $sample -d $depths -o ${plot_base_filename}-log_scale.pdf --log-scale-y
  plot_coverage.py --sample-name $sample -d $depths -v $vcf -o ${plot_base_filename}-with_variants.pdf
  plot_coverage.py --sample-name $sample -d $depths -v $vcf -o ${plot_base_filename}-with_variants-log_scale.pdf --log-scale-y
  plot_coverage.py --sample-name $sample --no-highlight -d $depths -o ${plot_base_filename}-no_low_cov_highlighting.pdf
  plot_coverage.py --sample-name $sample --no-highlight --log-scale-y -d $depths -o ${plot_base_filename}-log_scale-no_low_cov_highlighting.pdf
  plot_coverage.py --sample-name $sample --no-highlight -d $depths -v $vcf -o ${plot_base_filename}-with_variants-no_low_cov_highlighting.pdf
  plot_coverage.py --sample-name $sample --no-highlight --log-scale-y -d $depths -v $vcf -o ${plot_base_filename}-with_variants-log_scale-no_low_cov_highlighting.pdf
  """
}
