/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScflow.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.manifest ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
 * Create a channel for input read files
 */
if (params.manifest) { ch_manifest = file(params.manifest, checkIfExists: true) }
if (params.input) { ch_input = file(params.input, checkIfExists: true) }
if (params.ctd_path) { ch_ctd_path = file(params.ctd_path, checkIfExists: true) }
if (params.celltype_mappings) { ch_celltype_mappings = file(params.celltype_mappings, checkIfExists: false) }
if (params.ensembl_mappings) { ch_ensembl_mappings = file(params.ensembl_mappings, checkIfExists: false) }
if (params.reddim_genes_yml) { ch_reddim_genes_yml = file(params.reddim_genes_yml, checkIfExists: true) }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { SCFLOW_CHECKINPUTS         } from '../modules/local/process/scflow/checkinputs'       addParams( options: params.modules['scflow_checkinputs']     )
include { SCFLOW_QC                  } from '../modules/local/process/scflow/qc'                addParams( options: params.modules['scflow_qc']              )
include { SCFLOW_MERGEQCTABLES       } from '../modules/local/process/scflow/mergeqctables'     addParams( options: params.modules['scflow_mergeqctables']   )
include { SCFLOW_MERGE               } from '../modules/local/process/scflow/merge'             addParams( options: params.modules['scflow_merge']           )
include { SCFLOW_INTEGRATE           } from '../modules/local/process/scflow/integrate'         addParams( options: params.modules['scflow_integrate']       )
include { SCFLOW_REDUCEDIMS          } from '../modules/local/process/scflow/reducedims'        addParams( options: params.modules['scflow_reducedims']      )
include { SCFLOW_CLUSTER             } from '../modules/local/process/scflow/cluster'           addParams( options: params.modules['scflow_cluster']         )
include { SCFLOW_REPORTINTEGRATED    } from '../modules/local/process/scflow/reportintegrated'  addParams( options: params.modules['scflow_reportintegrated'])
include { SCFLOW_MAPCELLTYPES        } from '../modules/local/process/scflow/mapcelltypes'      addParams( options: params.modules['scflow_mapcelltypes']    )
include { SCFLOW_FINALIZE            } from '../modules/local/process/scflow/finalize'          addParams( options: params.modules['scflow_finalize']        )
include { SCFLOW_PLOTREDDIMGENES     } from '../modules/local/process/scflow/plotreddimgenes'   addParams( options: params.modules['scflow_plotreddimgenes'] )
include { SCFLOW_DGE                 } from '../modules/local/process/scflow/dge'               addParams( options: params.modules['scflow_dge']             )
include { SCFLOW_IPA                 } from '../modules/local/process/scflow/ipa'               addParams( options: params.modules['scflow_ipa']             )
include { SCFLOW_DIRICHLET           } from '../modules/local/process/scflow/dirichlet'         addParams( options: params.modules['scflow_dirichlet']       )
include { GET_SOFTWARE_VERSIONS      } from '../modules/local/get_software_versions'            addParams( options: [publish_files : ['tsv':'']]             )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

process UNZIP {
    tag "${key}"

    input:
    tuple val(key), path(mat_path)

    output:
    tuple val(key), path(mat_folder)

    script:
    """
    if [[ -d ${mat_path} ]]; then
        echo "${mat_path} is a directory"
        mv ${mat_path} mat_folder
    elif [[ -f ${mat_path} ]]; then
        echo "${mat_path} is a file"
        mkdir mat_folder && unzip ${mat_path} -d ./mat_folder
    else
        echo "${mat_path} is not valid"
        mv ${mat_path} mat_folder
        exit 1
    fi
    """
}

workflow SCFLOW {

    main:
    SCFLOW_CHECKINPUTS (
        ch_manifest,
        ch_input
    )

    UNZIP (
        SCFLOW_CHECKINPUTS.out.checked_manifest.splitCsv(
            header:['key', 'filepath'],
            skip: 1, sep: '\t'
            )
        .map { row -> tuple(row.key, row.filepath) },
        ch_input,
        ch_ensembl_mappings
    )

    SCFLOW_QC (
        UNZIP.out,
        ch_input,
        ch_ensembl_mappings
    )

    SCFLOW_MERGEQCTABLES (
        SCFLOW_QC.out.qc_summary.collect()
    )

    SCFLOW_MERGE (
        SCFLOW_QC.out.qc_sce.collect(),
        ch_ensembl_mappings
    )

    SCFLOW_INTEGRATE (
        SCFLOW_MERGE.out.merged_sce
    )

    SCFLOW_REDUCEDIMS (
        SCFLOW_INTEGRATE.out.integrated_sce
    )

    SCFLOW_CLUSTER (
        SCFLOW_REDUCEDIMS.out.reddim_sce
    )

    SCFLOW_REPORTINTEGRATED (
        SCFLOW_CLUSTER.out.clustered_sce
    )

    SCFLOW_MAPCELLTYPES (
        SCFLOW_CLUSTER.out.clustered_sce,
        ch_ctd_path
    )

    SCFLOW_FINALIZE (
        SCFLOW_MAPCELLTYPES.out.celltype_mapped_sce,
        ch_celltype_mappings
    )

    SCFLOW_DGE (
        SCFLOW_FINALIZE.out.final_sce,
        params.dge_de_method,
        SCFLOW_FINALIZE.out.celltypes.splitCsv (
            header:['celltype', 'n_cells'], skip: 1, sep: '\t'
        )
        .map { row -> tuple(row.celltype, row.n_cells) },
        ch_ensembl_mappings
    )

    SCFLOW_IPA (
        SCFLOW_DGE.out.de_table
    )

    SCFLOW_DIRICHLET (
        SCFLOW_FINALIZE.out.final_sce
    )

    SCFLOW_PLOTREDDIMGENES (
        SCFLOW_CLUSTER.out.clustered_sce,
        ch_reddim_genes_yml
    )

    GET_SOFTWARE_VERSIONS (
    )
}


/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
