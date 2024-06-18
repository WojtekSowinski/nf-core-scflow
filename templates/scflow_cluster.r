#!/usr/bin/env Rscript
# Reduce dimensions for a SCE
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

##  ............................................................................
##  Load packages                                                           ####
library(scFlow)
library(parallel)
library(SingleCellExperiment) # due to monocle3 missing namespace::
library(knitr) # due to missing knitr:: namespace in the integrate report

##  ............................................................................
##  Parse pipeline configuration

args <- {}
args$sce_path <- "!{sce_path}"
args$cluster_method <- "!{params.clust_cluster_method}"
args$reduction_method <- "!{params.clust_reduction_method}"
args$res <- !{params.clust_res}
args$k <- !{params.clust_k}
args$louvain_iter <- !{params.clust_louvain_iter}

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()

##  ............................................................................
##  Start                                                                   ####

sce <- read_sce(args$sce_path, read_metadata = TRUE)

sce <- cluster_sce(
  sce,
  cluster_method = args$cluster_method,
  reduction_method = args$reduction_method,
  res = args$res,
  k = args$k,
  louvain_iter = args$louvain_iter
)

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "clustered_sce"),
  write_metadata = TRUE
)

##  ............................................................................
##  Clean up                                                                ####

# Clear biomart cache
