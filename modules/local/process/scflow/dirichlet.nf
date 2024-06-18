/*
 * Dirichlet modeling of relative cell-type abundance
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_DIRICHLET {
    tag 'DIRICHLET'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //    container 'combiz/scflow-docker:0.6.1'

    input:
    path sce

    output:
    path 'dirichlet_report', emit: dirichlet_report

    shell:
    def software = getSoftwareName(task.process)

    template "scflow_dirichlet.r"

}
