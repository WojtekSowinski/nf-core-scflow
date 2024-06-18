#!/usr/bin/env Rscript
# Reduce dimensions for a SCE
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(parallel)
library(SingleCellExperiment) # due to monocle3 missing namespace::

##  ............................................................................
##  Parse pipeline configuration

args <- {}
args$sce_path <- "!{sce_path}"
args$input_reduced_dim <- "!{params.reddim_input_reduced_dim}"
args$reduction_methods <- "!{params.reddim_reduction_methods}"
args$vars_to_regress_out <- "!{params.reddim_vars_to_regress_out}"
args$pca_dims <- !{params.reddim_umap_pca_dims}
args$n_neighbors <- !{params.reddim_umap_n_neighbors}
args$n_components <- !{params.reddim_umap_n_components}
args$init <- "!{params.reddim_umap_init}"
args$metric <- "!{params.reddim_umap_metric}"
args$n_epochs <- !{params.reddim_umap_n_epochs}
args$learning_rate <- !{params.reddim_umap_learning_rate}
args$min_dist <- !{params.reddim_umap_min_dist}
args$spread <- !{params.reddim_umap_spread}
args$set_op_mix_ratio <- !{params.reddim_umap_set_op_mix_ratio}
args$local_connectivity <- !{params.reddim_umap_local_connectivity}
args$repulsion_strength <- !{params.reddim_umap_repulsion_strength}
args$negative_sample_rate <- !{params.reddim_umap_negative_sample_rate}
args$fast_sgd <- "!{params.reddim_umap_fast_sgd}"
args$dims <- !{params.reddim_tsne_dims}
args$initial_dims <- !{params.reddim_tsne_initial_dims}
args$perplexity <- !{params.reddim_tsne_perplexity}
args$theta <- !{params.reddim_tsne_theta}
args$stop_lying_iter <- !{params.reddim_tsne_stop_lying_iter}
args$mom_switch_iter <- !{params.reddim_tsne_mom_switch_iter}
args$max_iter <- !{params.reddim_tsne_max_iter}
args$pca_center <- "!{params.reddim_tsne_pca_center}"
args$pca_scale <- "!{params.reddim_tsne_pca_scale}"
args$normalize <- "!{params.reddim_tsne_normalize}"
args$momentum <- !{params.reddim_tsne_momentum}
args$final_momentum <- !{params.reddim_tsne_final_momentum}
args$eta <- !{params.reddim_tsne_eta}
args$exaggeration_factor <- !{params.reddim_tsne_exaggeration_factor}


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()
args$fast_sgd <- as.logical(args$fast_sgd)
args$input_reduced_dim <- strsplit(args$input_reduced_dim, ",")[[1]]
args$reduction_methods <- strsplit(args$reduction_methods, ",")[[1]]
args$vars_to_regress_out <- strsplit(args$vars_to_regress_out, ",")[[1]]
args <- purrr::map(args, function(x) {
  if (length(x) == 1) {
    if (toupper(x) == "TRUE") return(TRUE)
    if (toupper(x) == "FALSE") return(FALSE)
    if (toupper(x) == "NULL") return(NULL)
  }
  return(x)
})

##  ............................................................................
##  Start                                                                   ####

sce <- read_sce(args$sce_path, read_metadata = TRUE)

set.seed(42)
sce <- reduce_dims_sce(
  sce,
  input_reduced_dim = args$input_reduced_dim,
  reduction_methods = args$reduction_methods,
  vars_to_regress_out = args$vars_to_regress_out,
  pca_dims = args$pca_dims,
  n_neighbors = args$n_neighbors,
  n_components = args$n_components,
  init = args$init,
  metric = args$metric,
  n_epochs = args$n_epochs,
  learning_rate = args$learning_rate,
  min_dist = args$min_dist,
  spread = args$spread,
  set_op_mix_ratio = args$set_op_mix_ratio,
  local_connectivity = args$local_connectivity,
  repulsion_strength = args$repulsion_strength,
  negative_sample_rate = args$negative_sample_rate,
  fast_sgd = args$fast_sgd,
  dims = args$dims,
  initial_dims = args$initial_dims,
  perplexity = args$perplexity,
  theta = args$theta,
  stop_lying_iter = args$stop_lying_iter,
  mom_switch_iter = args$mom_switch_iter,
  max_iter = args$max_iter,
  pca_center = args$pca_center,
  pca_scale = args$pca_scale,
  pca_normalize = args$pca_normalize,
  momentum = args$momentum,
  final_momentum = args$final_momentum,
  eta = args$eta,
  exaggeration_factor = args$exaggeration_factor
)

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "reddim_sce"),
  write_metadata = TRUE
)

##  ............................................................................
##  Clean up                                                                ####

# Clear biomart cache
