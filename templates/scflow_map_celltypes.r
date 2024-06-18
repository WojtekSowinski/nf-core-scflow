#!/usr/bin/env Rscript
# Map celltypes for a SCE
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = !{task.cpus})
system("mkdir ctd_folder && unzip !{ctd_path} -d ./ctd_folder")

##  ............................................................................
##  Load packages                                                           ####
library(scFlow)
library(parallel)

##  ............................................................................
##  Parse pipeline configuration

args <- {}
args$sce_path <- "!{sce}"
args$ctd_folder <- "ctd_folder"
args$clusters_colname <- "!{params.cta_clusters_colname}"
args$cells_to_sample <- !{params.cta_cells_to_sample}
args$annotation_level <- !{params.cta_annotation_level}
args$species <- "!{params.species}"
args$reddimplot_pointsize <- !{params.reddimplot_pointsize}
args$reddimplot_alpha <- !{params.reddimplot_alpha}

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

options("scflow_species" = args$species)
options("scflow_reddimplot_pointsize" = args$reddimplot_pointsize)
options("scflow_reddimplot_alpha" = args$reddimplot_alpha)

##  ............................................................................
##  Start                                                                   ####

cat(print(tempdir()))

sce <- read_sce(args$sce_path)

sce <- map_celltypes_sce(
  sce,
  ctd_folder = args$ctd_folder,
  clusters_colname = args$clusters_colname,
  cells_to_sample = args$cells_to_sample,
  annotation_level = as.numeric(args$annotation_level),
  species = args$species
)

##  ............................................................................
##  Save Outputs                                                            ####

write_celltype_mappings(sce, folder_path = getwd())

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "celltype_mapped_sce")
)

scflow_version <- cat(as.character(utils::packageVersion("scFlow")))
cat("scFlow", scflow_version, file=paste0("scFlow_",scflow_version,".version.txt"))
