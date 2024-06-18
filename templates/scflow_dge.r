#!/usr/bin/env Rscript
# Perform differential gene expression on a SingleCellExperiment Object
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####


##  ............................................................................
##  Load packages                                                           ####
library(cli)
# Note: scFlow is loaded after the mc.cores option is defined/overriden below

##  ............................................................................
##  Parse pipeline configuration

args <- {}
args$ensembl_mappings <- "!{ensembl_mappings}"
args$de_method <- "!{de_method}"
args$celltype <- "!{celltype}"
args$sce <- "!{sce}"
args$mast_method <- "!{params.dge_mast_method}"
args$min_counts <- !{params.dge_min_counts}
args$min_cells_pc <- !{params.dge_min_cells_pc}
args$rescale_numerics <- "!{params.dge_rescale_numerics}"
args$force_run <- "!{params.dge_force_run}"
args$pseudobulk <- "!{params.dge_pseudobulk}"
args$celltype_var <- "!{params.dge_celltype_var}"
args$sample_var <- "!{params.dge_sample_var}"
args$dependent_var <- "!{params.dge_dependent_var}"
args$ref_class <- "!{params.dge_ref_class}"
args$confounding_vars <- "!{params.dge_confounding_vars}"
args$random_effects_var <- "!{params.dge_random_effects_var}"
args$padj_cutoff <- !{params.dge_padj_cutoff}
args$logFC_threshold <- !{params.dge_logFC_threshold}
args$species <- "!{params.species}"
args$max_cores <- "!{params.dge_max_cores}"

options("scflow_species" = args$species)

args$rescale_numerics <- as.logical(args$rescale_numerics)
args$pseudobulk <- as.logical(args$pseudobulk)
args$force_run <- as.logical(args$force_run)
if (tolower(args$random_effects_var) == "null") args$random_effects_var <- NULL

args$max_cores <- if (toupper(args$max_cores) == "NULL") {
  NULL
} else {
  as.numeric(as.character(args$max_cores))
}

args$confounding_vars <- strsplit(args$confounding_vars, ",")[[1]]

#   ____________________________________________________________________________
#   Delay Package Loading for Optional Max Cores Override

n_cores <- future::availableCores(methods = "mc.cores")

if (is.null(args$max_cores)) {
  options(mc.cores = n_cores)
} else {
  options(mc.cores = min(args$max_cores, n_cores))
}

cli::cli_alert(sprintf(
  "Using %s cores on system with %s available cores.",
  getOption("mc.cores"),
  n_cores
))


library(scFlow)

#   ____________________________________________________________________________
#   Start DE                                                                ####

write(sprintf(
  "##### Starting DE of %s cells with %s",
  args$celltype, args$demethod
), stdout())

sce <- read_sce(args$sce)

sce_subset <- sce[, sce$cluster_celltype == args$celltype]

if (args$pseudobulk) {
  pb_str <- "_pb"
  sce_subset <- pseudobulk_sce(
    sce_subset,
    keep_vars = c(
      args$dependent_var,
      args$confounding_vars,
      args$random_effects_var
    ),
    assay_name = "counts",
    celltype_var = args$celltype_var,
    sample_var = args$sample_var
  )
} else {
  pb_str <- ""
}

de_results <- perform_de(
  sce_subset,
  de_method = args$de_method,
  min_counts = args$min_counts,
  min_cells_pc = args$min_cells_pc,
  rescale_numerics = args$rescale_numerics,
  dependent_var = args$dependent_var,
  ref_class = args$ref_class,
  confounding_vars = args$confounding_vars,
  random_effects_var = args$random_effects_var,
  mast_method = args$mast_method,
  force_run = args$force_run,
  ensembl_mapping_file = args$ensembl_mappings,
  species = getOption("scflow_species")
)

file_name <- paste0(
  args$celltype, "_",
  args$de_method, pb_str, "_"
)

for (result in names(de_results)) {
  if (dim(de_results[[result]])[[1]] > 0) {
    write.table(de_results[[result]],
      file = file.path(
        getwd(),
        paste0(file_name, result, "_DE.tsv")
      ),
      quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
    )

    report_de(de_results[[result]],
      logFC_threshold = args$logFC_threshold,
      padj_cutoff = args$padj_cutoff,
      n_label = args$n_label,
      report_folder_path = file.path(getwd()),
      report_file = paste0(file_name, result, "_scflow_de_report")
    )

    print("report generated")

    p <- scFlow::volcano_plot(
      dt = de_results[[result]],
      logFC_threshold = args$logFC_threshold,
      padj_cutoff = args$padj_cutoff,
      n_label = args$n_label
    )

    ggplot2::ggsave(
      filename = file.path(
        getwd(),
        paste0(file_name, result, "_volcano_plot.png")
      ),
      plot = p,
      width = 7, height = 5, units = "in", dpi = 600
    )

    print("Volcano plot generated")
  } else {
    print(sprintf("No DE genes found for %s", result))
  }
}
