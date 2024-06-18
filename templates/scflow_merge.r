#!/usr/bin/env Rscript
# Merge multiple SingleCellExperiments
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

##  ............................................................................
##  Load packages                                                           ####
library(scFlow)
library(SingleCellExperiment) # due to monocle3 missing namespace::

##  ............................................................................
##  Parse configuration

args <- {}
args$unique_id_var <- "!{params.qc_key_colname}"
args$plot_vars <- "!{params.merge_plot_vars}"
args$facet_vars <- "!{params.merge_facet_vars}"
args$outlier_vars <- "!{params.merge_outlier_vars}"
args$species <- "!{params.species}"
args$sce_paths <- "!{sce_paths}"
args$ensembl_mappings <- "!{ensembl_mappings}"

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

options("scflow_species" = args$species)

args$sce_paths <- strsplit(args$sce_paths, ",")[[1]]
args$facet_vars <- strsplit(args$facet_vars, ",")[[1]]
args$outlier_vars <- strsplit(args$outlier_vars, ",")[[1]]
args$plot_vars <- strsplit(args$plot_vars, ",")[[1]]

args <- purrr::map(args, function(x) {
  if (length(x) == 1) {
    if (toupper(x) == "TRUE") return(TRUE)
    if (toupper(x) == "FALSE") return(FALSE)
    if (toupper(x) == "NULL") return(NULL)
  }
  return(x)
})


##  ............................................................................
##  Start Merge                                                             ####

print(sprintf(
  "Reading %sx SingleCellExperiment's",
  length(args$sce_paths))
)

sce_l <- lapply(args$sce_paths, read_sce)

sce <- merge_sce(
  sce_l,
  ensembl_mapping_file = args$ensembl_mappings
)

dir.create(file.path(getwd(), "merged_report"))

sce <- annotate_merged_sce(
  sce,
  plot_vars = args$plot_vars,
  unique_id_var = args$unique_id_var,
  facet_vars = args$facet_vars,
  outlier_vars = args$outlier_vars
)

report_merged_sce(
  sce = sce,
  report_folder_path = file.path(getwd(), "merged_report"),
  report_file = paste0(args$unique_id_var, "_scflow_merged_report")
)

##  ............................................................................
##  Save Outputs                                                            ####

dir.create(file.path(getwd(), "merge_plots"))
# Save merged plots (images)
for (rd_name in setdiff(names(sce@metadata$pseudobulk_rd_plots), "UMAP3D")) {
  png(file.path(getwd(), "merge_plots",
                paste0(args$unique_id_var, "_", rd_name, ".png")),
      width = 247, height = 170, units = "mm", res = 600)
  print(sce@metadata$pseudobulk_rd_plots[[rd_name]])
  dev.off()
}

dir.create(file.path(getwd(), "pb_plots"))
# Save pb plots (images)
for (rd_name in names(sce@metadata$pseudobulk_plots)) {
  png(file.path(getwd(), "pb_plots",
                paste0(args$unique_id_var, "_", rd_name, ".png")),
      width = 247, height = 170, units = "mm", res = 600)
  print(sce@metadata$pseudobulk_rd_plots[[rd_name]])
  dev.off()
}

dir.create(file.path(getwd(), "merge_summary_plots"))
# save multi-sample summary plots
for (plot_var in names(sce@metadata$merged_plots)) {
  for (plot in names(sce@metadata$merged_plots[[plot_var]])) {
    if (plot_var == plot) {
      plot_caption <- sprintf(
        "%s_by_%s",
        plot_var, sce@metadata$merge_qc_params$unique_id_var)
    } else {
      plot_caption <- sprintf(
        "%s_by_%s",
        plot_var, strsplit(plot, "_vs_")[[1]][[2]])
    }
    png(file.path(getwd(), "merge_summary_plots",
                  paste0(args$unique_id_var, "_", plot_caption, ".png")),
        width = 247, height = 170, units = "mm", res = 600)
    print(sce@metadata$merged_plots[[plot_var]][[plot]])
    dev.off()
  }
}

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "merged_sce")
)


##  ............................................................................
##  Clean up                                                                ####
