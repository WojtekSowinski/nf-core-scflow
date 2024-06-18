/*
 * Impacted pathway analysis of differentially expressed genes
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_IPA {
    tag "${de_table_basename}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //    container 'combiz/scflow-docker:0.6.1'

    input:
    path de_table

    output:
    path '*_ipa'      , emit: ipa_results , optional: true, type: 'dir'
    path '*.html'   , emit: ipa_report  , optional: true

    shell:
    de_table_basename = "${de_table.baseName}"
    def software = getSoftwareName(task.process)

    template "scflow_ipa.r"

}
