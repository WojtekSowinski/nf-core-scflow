#!/usr/bin/env Rscript
# Model relative celltype abundance with a Dirichlet model
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = !{taks.cpus})

##  ............................................................................
##  Load packages                                                           ####
library(scFlow)
library(cli)

##  ............................................................................
##  Parse pipeline configuration

args <- {}
args$sce_path <- "!{sce}"
args$unique_id_var <- "!{params.dirich_unique_id_var}"
args$celltype_var <- "!{params.dirich_celltype_var}"
args$dependent_var <- "!{params.dirich_dependent_var}"
args$ref_class <- "!{params.dirich_ref_class}"
args$var_order <- "!{params.dirich_var_order}"
args$confounding_vars <- "!{params.dirich_confounding_vars}"

if (tolower(args$var_order) == "null") {
  args$var_order <- NULL
  } else {
    args$var_order <- strsplit(args$var_order, ",")[[1]]
  }

if (tolower(args$confounding_vars) == "null") {
  args$confounding_vars <- NULL
  } else {
    args$confounding_vars <- strsplit(args$confounding_vars, ",")[[1]]
  }

#   ____________________________________________________________________________
#   Start                                                                   ####

sce <- read_sce(args$sce_path)

results <- model_celltype_freqs(
  sce,
  unique_id_var = args$unique_id_var,
  celltype_var = args$celltype_var,
  dependent_var = args$dependent_var,
  ref_class = args$ref_class,
  var_order = args$var_order,
  confounding_vars = args$confounding_vars
)

## ............................................................................
## Save Outputs ####

new_dirs <- c(
  "dirichlet_report"
)

#make dirs
purrr::walk(new_dirs, ~ dir.create(file.path(getwd(), .)))

report_celltype_model(
  results,
  report_folder_path = file.path(
    getwd(),
    "dirichlet_report"
  ),
  report_file = paste0(
    args$celltype_var,
    args$dependent_var,
    "dirichlet_report",
    sep = "_")
)

scflow_version <- cat(as.character(utils::packageVersion("scFlow")))
cat("scFlow", scflow_version, file=paste0("scFlow_",scflow_version,".version.txt"))
