#!/usr/bin/env Rscript
# Reduce dimensions for a SCE
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)

##  ............................................................................
##  Parse command-line arguments                                            ####

args <- {}
args$sce_path <- "!{sce}"
args$reddim_genes_yml <- "!{reddim_genes_yml}"
args$reduction_methods <- "!{params.plotreddim_reduction_methods}"
args$reddimplot_pointsize <- !{params.reddimplot_pointsize}
args$reddimplot_alpha <- !{params.reddimplot_alpha}

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args$reduction_methods <- strsplit(args$reduction_methods, ",")[[1]]
options("scflow_reddimplot_pointsize" = args$reddimplot_pointsize)
options("scflow_reddimplot_alpha" = args$reddimplot_alpha)

##  ............................................................................
##  Start                                                                   ####

sce <- read_sce(args$sce_path)

gene_l <- yaml::read_yaml(args$reddim_genes_yml)

valid_reddims <- SingleCellExperiment::reducedDimNames(sce)

assertthat::assert_that(
  all(args$reduction_methods %in% valid_reddims),
  msg = sprintf("Valid reddims are: %s", paste0(valid_reddims, collapse = ",")))

for (reddim in args$reduction_method) {
  for (l_name in names(gene_l)) {
    folder_path <- file.path(getwd(), "reddim_gene_plots", reddim, l_name)
    R.utils::mkdirs(folder_path)
    for (gene in gene_l[[l_name]]) {
      print(gene)
      if (gene %in% SummarizedExperiment::rowData(sce)$gene) {
        p <- plot_reduced_dim_gene(
          sce,
          reduced_dim = reddim,
          gene = gene
        )
        png(file.path(folder_path, paste0(gene, ".png")),
            width = 170, height = 170, units = "mm", res = 600)
        print(p)
        dev.off()
      } else {
        warning(print(sprintf("Gene %s not found.", gene)))
      }
    }
  }
}

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment


##  ............................................................................
##  Clean up                                                                ####

# Clear biomart cache
