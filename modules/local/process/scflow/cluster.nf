/*
 * CLUSTERING
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_CLUSTER {
    tag 'MERGED'
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //container 'combiz/scflow-docker:0.6.1'

    input:
        path sce

    output:
        path 'clustered_sce/'       , emit: clustered_sce, type: 'dir'

    shell:
    template "scflow_cluster.r"

}
