/*
 * Generate integration HTML report
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_REPORTINTEGRATED {
    tag 'MERGED'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //    container 'combiz/scflow-docker:0.6.1'

    input:
    path( sce )

    output:
    path 'integration_report', emit: integration_report, type: 'dir'

    shell:
    def software = getSoftwareName(task.process)

    template "scflow_report_integrated.r"

}
