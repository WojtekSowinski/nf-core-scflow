/*
 * Annotate cluster celltypes
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_MAPCELLTYPES {
    tag 'MERGED'
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //    container 'combiz/scflow-docker:0.6.1'

    input:
    path sce
    path ctd_path

    output:
    path 'celltype_mapped_sce/' , emit: celltype_mapped_sce, type: 'dir'
    path 'celltype_mappings.tsv', emit: celltype_mappings

    shell:
    def software = getSoftwareName(task.process)

    template "scflow_map_celltypes.r"

}
