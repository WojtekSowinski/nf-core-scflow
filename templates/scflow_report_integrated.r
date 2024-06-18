#!/usr/bin/env Rscript
# Reduce dimensions for a SCE
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = !{task.cpus});

##  ............................................................................
##  Load packages                                                           ####
library(scFlow)
library(parallel)
library(SingleCellExperiment) # due to monocle3 missing namespace::
library(knitr) # due to missing knitr:: namespace in the integrate report

##  ............................................................................
##  Parse pipeline configuration

args <- {}
args$sce_path <- "!{sce}"
args$categorical_covariates <- "!{params.integ_categorical_covariates}"
args$input_reduced_dim <- "!{params.integ_input_reduced_dim}"
args$reddimplot_pointsize <- !{params.reddimplot_pointsize}
args$reddimplot_alpha <- !{params.reddimplot_alpha}

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args$categorical_covariates <- strsplit(args$categorical_covariates, ",")[[1]]

options("scflow_reddimplot_pointsize" = args$reddimplot_pointsize)
options("scflow_reddimplot_alpha" = args$reddimplot_alpha)


##  ............................................................................
##  Start                                                                   ####

sce <- read_sce(args$sce_path, read_metadata = TRUE)

sce <- annotate_integrated_sce(
  sce,
  categorical_covariates = args$categorical_covariates,
  input_reduced_dim = args$input_reduced_dim
)

##  ............................................................................
##  Save Outputs                                                            ####

dir.create(file.path(getwd(), "integration_report"))

report_integrated_sce(
  sce = sce,
  report_folder_path = file.path(getwd(), "integration_report"),
  report_file = "integrate_report_scflow"
)

scflow_version <- cat(as.character(utils::packageVersion("scFlow")))
cat("scFlow", scflow_version, file=paste0("scFlow_",scflow_version,".version.txt"))
