#!/usr/bin/env Rscript
# Perform quality-control on a feature-barcode matrix with scflow
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = !{task.cpus})

##  ............................................................................
##  Load packages                                                           ####
library(scFlow)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- {}
args$key_colname <- "!{params.qc_key_colname}"
args$factor_vars <- "!{params.qc_factor_vars}"
args$min_library_size <- !{params.qc_min_library_size}
args$max_library_size <- "!{params.qc_max_library_size}"
args$min_features <- !{params.qc_min_features}
args$max_features <- "!{params.qc_max_features}"
args$max_mito <- "!{params.qc_max_mito}"
args$min_ribo <- !{params.qc_min_ribo}
args$max_ribo <- !{params.qc_max_ribo}
args$min_counts <- !{params.qc_min_counts}
args$min_cells <- !{params.qc_min_cells}
args$drop_unmapped <- "!{params.qc_drop_unmapped}"
args$drop_mito <- "!{params.qc_drop_mito}"
args$drop_ribo <- "!{params.qc_drop_ribo}"
args$nmads <- !{params.qc_nmads}
args$find_singlets <- "!{params.mult_find_singlets}"
args$singlets_method <- "!{params.mult_singlets_method}"
args$vars_to_regress_out <- "!{params.mult_vars_to_regress_out}"
args$pca_dims <- !{params.mult_pca_dims}
args$var_features <- !{params.mult_var_features}
args$doublet_rate <- !{params.mult_doublet_rate}
args$dpk <- "!{params.mult_dpk}"
args$pK <- "!{params.mult_pK}"
args$find_cells <- "!{params.amb_find_cells}"
args$lower <- !{params.amb_lower}
args$retain <- "!{params.amb_retain}"
args$alpha_cutoff <- !{params.amb_alpha_cutoff}
args$niters <- !{params.amb_niters}
args$expect_cells <- "!{params.amb_expect_cells}"
args$species <- "!{params.species}"


options("scflow_species" = args$species)

args[startsWith(names(args), "drop_")] <-
  as.logical(args[startsWith(names(args), "drop_")])
args$max_library_size <- ifelse(
  args$max_library_size == "adaptive",
  args$max_library_size,
  as.numeric(as.character(args$max_library_size))
)
args$max_features <- ifelse(
  args$max_features == "adaptive",
  args$max_features,
  as.numeric(as.character(args$max_features))
)
args$max_mito <- ifelse(
  args$max_mito == "adaptive",
  args$max_mito,
  as.numeric(as.character(args$max_mito))
)
args$pK <- if (toupper(args$pK) == "NULL") NULL else {
  as.numeric(as.character(args$pK))
}
args$dpk <- if (toupper(args$dpk) == "NULL") NULL else {
  as.numeric(as.character(args$dpk))
}

if (toupper(args$retain) == "NULL") {
  args$retain <- NULL
} else if (toupper(args$retain) == "AUTO") {
  args$retain <- "auto"
} else {
  args$retain <- as.numeric(as.character(args$retain))
}

args$find_singlets <- as.logical(args$find_singlets)
args$vars_to_regress_out <- strsplit(args$vars_to_regress_out, ",")[[1]]
args <- purrr::map(args, function(x) {
  if (length(x) == 1) {
    if (toupper(x) == "TRUE") return(TRUE)
    if (toupper(x) == "FALSE") return(FALSE)
    if (toupper(x) == "NULL") return(NULL)
  }
  return(x)
})

if (!is.null(args$factor_vars)) {
  args$factor_vars <- strsplit(args$factor_vars, ",")[[1]]
  col_classes <- rep("factor", length(args$factor_vars))
  names(col_classes) <- args$factor_vars
} else {
  col_classes <- NA
}



##  ............................................................................
##  Start QC                                                                ####

cli::boxx("Analysing: !{key}", float = "center")

mat <- scFlow::read_sparse_matrix("!{mat_path}")

metadata <- read_metadata(
  unique_key = "!{key}",
  key_colname = args$key_colname,
  samplesheet_path = "!{input}",
  col_classes = col_classes
)

sce <- generate_sce(mat, metadata)

rm(mat, metadata)

if (args$find_cells) {
  sce <- find_cells(
    sce,
    lower = args$lower,
    retain = args$retain,
    alpha_cutoff = args$alpha_cutoff,
    niters = args$niters
  )
}

sce <- annotate_sce(
  sce = sce,
  min_library_size = args$min_library_size,
  max_library_size = args$max_library_size,
  min_features = args$min_features,
  max_features = args$max_features,
  max_mito = args$max_mito,
  min_ribo = args$min_ribo,
  max_ribo = args$max_ribo,
  min_counts = args$min_counts,
  min_cells = args$min_cells,
  drop_unmapped = args$drop_unmapped,
  drop_mito = args$drop_mito,
  drop_ribo = args$drop_ribo,
  nmads = args$nmads,
  annotate_genes = TRUE,
  annotate_cells = TRUE,
  ensembl_mapping_file = "!{ensembl_mappings}",
  species = args$species
)

sce <- filter_sce(
  sce,
  filter_genes = TRUE,
  filter_cells = TRUE
)

if (args$find_singlets) {
  sce <- find_singlets(
    sce = sce,
    singlet_find_method = args$singlets_method,
    vars_to_regress_out = args$vars_to_regress_out,
    pca_dims = args$pca_dims,
    var_features = args$var_features,
    doublet_rate = args$doublet_rate,
    dpk = args$dpk,
    pK = args$pK,
    num.cores = future::availableCores()
  )
  sce <- filter_sce(
    sce,
    filter_genes = TRUE,
    filter_cells = TRUE
  )
}


sce <- sce[ , sce$total_counts >= args$min_library_size]
sce <- sce[ , sce$total_features_by_counts >= args$min_features]


dir.create(file.path(getwd(), "qc_report"))

report_qc_sce(
  sce = sce,
  #report_folder_path = file.path(getwd(), "qc_report"),
  report_folder_path = file.path(getwd()),
  report_file = "!{key}_scflow_qc_report"
)

print("Analysis complete, saving outputs..")

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "!{key}_sce")
)

new_dirs <- c(
  "qc_plot_data",
  "qc_summary",
  "qc_plots")

#make dirs
purrr::walk(new_dirs, ~ dir.create(file.path(getwd(), .)))

# Save QC plots (tables)
for (df in names(sce@metadata$qc_plot_data)) {
  write.table(
    sce@metadata$qc_plot_data[[df]],
    file.path(getwd(), "qc_plot_data",
              paste0("!{key}", "_", df, ".tsv")),
    sep = "\t",
    col.names = TRUE, row.names = FALSE)
}

# Save QC summary table
write.table(
  cbind(sce@metadata$metadata, sce@metadata$qc_summary),
  file.path(getwd(), "qc_summary", "!{key}_qc_summary.tsv"),
  sep = "\t",
  col.names = TRUE, row.names = FALSE)

# Save QC plots (images)
for (pname in names(sce@metadata$qc_plots)) {
  png(file.path(getwd(), "qc_plots",
                paste0("!{key}", "_", pname, ".png")),
      width = 247, height = 170, units = "mm", res = 600)
  print(sce@metadata$qc_plots[[pname]])
  dev.off()
}

# Save doublet finder plots, square
for (pname in names(sce@metadata$qc_plots$doublet_finder)) {
  png(file.path(getwd(), "qc_plots",
                paste0("!{key}", "_", pname, "_doublet_finder.png")),
      width = 170, height = 170, units = "mm", res = 600)
  print(sce@metadata$qc_plots$doublet_finder[[pname]])
  dev.off()
}

system("mkdir sce")
system("mv !{key}_sce sce/")
scflow_version <- cat(as.character(utils::packageVersion("scFlow")))
cat("scFlow", scflow_version, file=paste0("scFlow_",scflow_version,".version.txt"))
