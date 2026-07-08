#!/usr/bin/env Rscript
# =============================================================================
# RIMA cross-atlas analysis (Phase 2)
#
# Maps brain-organoid cell-state neighbourhoods onto fetal-brain neighbourhoods
# to quantify transcriptomic maturation. Consumes the scVI-harmonized .h5ad
# files produced by the Python harmonization stage (run_harmonization.py) and
# writes match tables, conserved-gene tables, and plots to the --output dir.
#
# Usage:
#   Rscript run_rima.R --input <dir with harmonized h5ad files> --output <results dir>
#
# The stage functions live in sibling files (data_io.R, neighborhoods.R,
# matching.R, plots.R, cell_extraction.R, maturation.R, dtw_alignment.R) which
# are sourced below.
# =============================================================================

# --- Imports ---
library(RIMA)
library(miloR)
library(SingleCellExperiment)
library(anndataR)
library(reticulate)
library(ggplot2)
reticulate::py_config()

# --- Locate and source the pipeline modules (co-located with this script) ---
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}
script_dir <- get_script_dir()
source(file.path(script_dir, "data_io.R"))
source(file.path(script_dir, "neighborhoods.R"))
source(file.path(script_dir, "matching.R"))
source(file.path(script_dir, "plots.R"))
source(file.path(script_dir, "cell_extraction.R"))
source(file.path(script_dir, "maturation.R"))
source(file.path(script_dir, "dtw_alignment.R"))

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[which(args == "--input") + 1]
output_dir <- args[which(args == "--output") + 1]
if (length(input_dir) == 0 || length(output_dir) == 0) {
  stop("Missing --input or --output arguments")
}

# --- Step 0: Load harmonized datasets as SingleCellExperiments ---
datasets <- load_datasets(input_dir)
sce_fetal <- datasets$fetal
sce_organoid <- datasets$organoid

# --- Step 1: Build Milo neighbourhoods on the scVI latent space ---
print("Defining neighborhoods...")
mi_fetal <- define_neighbourhoods(sce_fetal, prop_seeds = 0.005, knn = 30)
mi_organoid <- define_neighbourhoods(sce_organoid, prop_seeds = 0.005, knn = 30)

print("Preprocessing Milo objects...")
milos <- preprocess_milos(mi_fetal, mi_organoid)
inspect_milos(milos)

# --- Step 2: Neighbourhood similarities ---
dt_sims <- calculate_similarities(milos, method = "spearman")

# Patch missing cell types before scrambling (returns the modified list).
milos <- patch_missing_celltypes(milos, col = "unified_celltype")

# --- Step 3: Significance -> matching -> conserved genes ---
dt_sims_sig <- calculate_significance(milos, dt_sims)
dt_match <- match_significant_nhoods(dt_sims_sig)
dt_cope <- calculate_conserved_genes(milos, dt_match, genes = NULL)

# --- Step 4: Save result tables ---
print("Saving tables to output directory...")
write.csv(dt_match, file.path(output_dir, "significant_matches.csv"), row.names = FALSE)
write.csv(dt_cope, file.path(output_dir, "top_conserved_genes.csv"), row.names = FALSE)

# --- Step 5: Save analysis plots (multi-page PDF) ---
save_analysis_plots(
  milos, dt_match, dt_cope,
  sce_fetal, sce_organoid, mi_fetal, mi_organoid,
  file.path(output_dir, "RIMA_Analysis_Plots.pdf")
)

# --- Step 6: Extract individual cells from the top match ---
write_match_summary(1, dt_match, mi_fetal, mi_organoid, sce_fetal, sce_organoid, output_dir)

# --- Step 7: Maturation age transfer + cell-type-specific plot ---
fetal_nhood_ages <- compute_fetal_nhood_ages(mi_fetal, sce_fetal)
mapped_ages <- build_mapped_ages(dt_match, fetal_nhood_ages, mi_organoid, sce_organoid)
maturation_plot <- plot_maturation_by_celltype(mapped_ages)
ggsave(file.path(output_dir, "Organoid_Maturation_by_CellType.pdf"),
       plot = maturation_plot, width = 14, height = 10)

# --- Step 8: DTW temporal alignment (one PDF per cell type) ---
run_dtw_alignment(dt_sims_sig, fetal_nhood_ages, mi_organoid, sce_organoid, output_dir)

print("Pipeline complete! Lineage-specific DTW plots saved successfully.")
