/*
 * Generate final SCE with optionally revised cell-types
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_FINALIZE {
    tag 'MERGED'
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //    container 'combiz/scflow-docker:0.6.1'

    input:
    path sce
    path celltype_mappings

    output:
    path 'final_sce'                , emit: final_sce, type: 'dir'
    path 'celltypes.tsv'            , emit: celltypes
    path 'celltype_metrics_report'  , emit: celltype_metrics_report, type: 'dir'
    path 'celltype_marker_tables'   , emit: celltype_marker_tables, type: 'dir'
    path 'celltype_marker_plots'    , emit: celltype_marker_plots, type: 'dir'

    shell:
    def software = getSoftwareName(task.process)
    def ctm = celltype_mappings.simpleName != 'NO_FILE' ? "$celltype_mappings" : 'nofile'

    template "scflow_finalize_sce.r"
}
