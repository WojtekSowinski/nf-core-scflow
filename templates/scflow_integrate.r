#!/usr/bin/env Rscript
# Integrate multiple single cell datasets (samples)
# Mahdi Moradi Marjaneh

# ____________________________________________________________________________
# Initialization ####

options(mc.cores = max(2, future::availableCores(methods = "mc.cores")))

## ............................................................................
## Load packages ####
library(scFlow)
library(parallel)

## ............................................................................
## Parse command-line arguments ####

args <- {}
args$sce_path <- "!{sce_path}"
args$method <- "!{params.integ_method}"
args$k <- !{params.integ_k}
args$unique_id_var <- "!{params.integ_unique_id_var}"
args$take_gene_union <- "!{params.integ_take_gene_union}"
args$remove_missing <- "!{params.integ_remove_missing}"
args$num_genes <- !{params.integ_num_genes}
args$combine <- "!{params.integ_combine}"
args$capitalize <- "!{params.integ_capitalize}"
args$use_cols <- "!{params.integ_use_cols}"
args$lambda <- !{params.integ_lambda}
args$thresh <- !{params.integ_thresh}
args$max_iters <- !{params.integ_max_iters}
args$nrep <- !{params.integ_nrep}
args$rand_seed <- !{params.integ_rand_seed}
args$quantiles <- !{params.integ_quantiles}
args$ref_dataset <- "!{params.integ_ref_dataset}"
args$min_cells <- !{params.integ_min_cells}
args$knn_k <- !{params.integ_knn_k}
args$center <- "!{params.integ_center}"
args$resolution <- !{params.integ_resolution}


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args ####

args <- purrr::map(args, function(x) {
  if (length(x) == 1) {
    if (toupper(x) == "TRUE") {
      return(TRUE)
    }
    if (toupper(x) == "FALSE") {
      return(FALSE)
    }
    if (toupper(x) == "NULL") {
      return(NULL)
    }
  }
  return(x)
})


## ............................................................................
## Integrate sce ####

sce <- read_sce(args$sce_path)

sce <- integrate_sce(
  sce,
  method = args$method,
  unique_id_var = args$unique_id_var,
  take_gene_union = args$take_gene_union,
  remove.missing = args$remove_missing,
  num_genes = args$num_genes,
  combine = args$combine,
  capitalize = args$capitalize,
  use_cols = args$use_cols,
  num_cores = future::availableCores(methods = "mc.cores"),
  k = args$k,
  lambda = args$lambda,
  thresh = args$thresh,
  max_iters = args$max_iters,
  nrep = args$nrep,
  H_init = NULL,
  W_init = NULL,
  V_init = NULL,
  rand_seed = args$rand_seed,
  knn_k = args$knn_k,
  ref_dataset = args$ref_dataset,
  min_cells = args$min_cells,
  quantiles = args$quantiles,
  resolution = args$resolution,
  center = args$center,
  print_obj = FALSE
)


## ............................................................................
## Save Outputs ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "integrated_sce"),
  write_metadata = TRUE
)
