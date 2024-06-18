/*
 * Perform dimensionality reduction
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_REDUCEDIMS {
    tag 'MERGED'
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //    container 'combiz/scflow-docker:0.6.1'

    input:
    path sce

    output:
    path 'reddim_sce/', emit: reddim_sce

    shell:
    def software = getSoftwareName(task.process)

    template "scflow_reduce_dims.r"

}
