/*
 * Run differential gene expression analysis
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_DGE {
    tag "${celltype} (${n_cells_str} cells) | ${de_method}"
    label 'process_high'
    errorStrategy 'ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"${celltype}_${de_method}") }

    //    container 'combiz/scflow-docker:0.6.1'

    input:
    path sce
    each de_method
    each ct_tuple
    path ensembl_mappings

    output:
    path '*.tsv'              , emit: de_table      , optional: true
    path '*.html'             , emit: de_report     , optional: true
    path '*_volcano_plot.png' , emit: de_plot       , optional: true
    //path 'de_plot_data/'               , emit: de_plot_data , optional: true

    shell:
    celltype     = ct_tuple[0]
    n_cells      = ct_tuple[1].toInteger()
    n_cells_str  = (Math.round(n_cells * 100) / 100000).round(1).toString() + 'k'
    def software = getSoftwareName(task.process)

    template "scflow_dge.r"
}
