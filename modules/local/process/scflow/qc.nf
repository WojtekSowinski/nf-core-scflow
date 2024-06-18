/*
 * Single Sample QC
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_QC {
    tag "${key}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"${key}") }

    input:
    tuple val(key), path(mat_path)
    path input
    path ensembl_mappings

    output:
    path '*.html'               , emit: qc_report
    path 'qc_plot_data'         , emit: qc_plot_data, type: 'dir'
    path 'qc_summary/*.tsv'     , emit: qc_summary
    path 'qc_plots'             , emit: qc_plots, type: 'dir'
    path 'sce/*_sce'            , emit: qc_sce, type: 'dir'

    shell:
    template "scflow_qc.r"
}
