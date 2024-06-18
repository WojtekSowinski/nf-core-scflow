#!/usr/bin/env Rscript
# Finalize SCE with manually revised celltypes
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

##  ............................................................................
##  Parse pipeline configuration

args <- {}
args$sce_path <- "!{sce_path}"
args$celltype_mappings <- "!{celltype_mappings}"
args$clusters_colname <- "!{params.cta_clusters_colname}"
args$celltype_var <- "!{params.cta_celltype_var}"
args$unique_id_var <- "!{params.cta_unique_id_var}"
args$facet_vars <- "!{params.cta_facet_vars}"
args$input_reduced_dim <- "!{params.clust_reduction_method}"
args$metric_vars <- "!{params.cta_metric_vars}"
args$top_n <- !{params.cta_top_n}
args$reddimplot_pointsize <- !{params.reddimplot_pointsize}
args$reddimplot_alpha <- !{params.reddimplot_alpha}
args$max_cores <- "!{params.max_cores}"

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args$facet_vars <- strsplit(args$facet_vars, ",")[[1]]
args$metric_vars <- strsplit(args$metric_vars, ",")[[1]]

options("scflow_reddimplot_pointsize" = args$reddimplot_pointsize)
options("scflow_reddimplot_alpha" = args$reddimplot_alpha)

args$max_cores <- if (toupper(args$max_cores) == "NULL") {
  NULL
} else {
  as.numeric(as.character(args$max_cores))
}

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
library(magrittr)
library(SingleCellExperiment)

##  ............................................................................
##  Start                                                                   ####

sce <- read_sce(args$sce_path)

if (file.exists(args$celltype_mappings)) {
  celltype_mappings <- read_celltype_mappings(args$celltype_mappings)
  sce <- map_custom_celltypes(
    sce,
    mappings = celltype_mappings,
    clusters_colname = args$clusters_colname
  )
} else {
  print("Revised cell-type mappings not provided, using auto-annotations")
}

sce <- annotate_celltype_metrics(
  sce,
  cluster_var = args$clusters_colname,
  celltype_var = args$celltype_var,
  unique_id_var = args$unique_id_var,
  facet_vars = args$facet_vars,
  input_reduced_dim = args$input_reduced_dim,
  metric_vars = args$metric_vars,
  top_n = args$top_n
)

dir.create(file.path(getwd(), "celltype_metrics_report"))

report_celltype_metrics(
  sce = sce,
  report_folder_path = file.path(getwd(), "celltype_metrics_report"),
  report_file = "scflow_celltype_metrics_report"
)

##  ............................................................................
##  Save Outputs                                                            ####

### Save cell-types/n_cells for NextFlow tags
celltypes <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
  dplyr::count(cluster_celltype)
colnames(celltypes) <- c("celltype", "n_cells")

write.table(
  data.frame(celltypes),
  file = "celltypes.tsv",
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

### Save Marker Gene Plots
folder_path <- file.path(getwd(), "celltype_marker_plots")
dir.create(folder_path)

for (group in names(sce@metadata$markers)) {
  pwidth <- max(
    10,
    length(unique(sce@metadata$markers[[group]]$marker_plot$data$Group))
  )
  pheight <- length(unique(sce@metadata$markers[[group]]$marker_plot$data$Gene))

  p <- sce@metadata$markers[[group]]$marker_plot

  plot_file_name <- paste0(group, "_markers")

  # save PNG
  png(file.path(folder_path, paste0(plot_file_name, ".png")),
    width = pwidth * 12, height = pheight * 5, units = "mm", res = 600
  )
  print(p)
  dev.off()

  # save PDF
  ggsave(
    file.path(folder_path, paste0(group, ".pdf")),
    p,
    width = pwidth * 12,
    height = pheight * 5,
    units = "mm",
    scale = 1
  )
}

### Save Marker Gene Tables
folder_path <- file.path(getwd(), "celltype_marker_tables")
dir.create(folder_path)
for (group in names(sce@metadata$markers)) {
  marker_test_file_name <- paste0(group, "_markers_test.tsv")
  top_markers_file_name <- paste0(group, "_top_markers.tsv")

  write.table(
    sce@metadata$markers[[group]]$marker_test_res,
    file = file.path(folder_path, marker_test_file_name),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )

  write.table(
    sce@metadata$markers[[group]]$top_specific_markers,
    file = file.path(folder_path, top_markers_file_name),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
}


### Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "final_sce")
)

