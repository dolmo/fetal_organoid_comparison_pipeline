# Imports 
library(RIMA)
library(miloR)
library(SingleCellExperiment)
library(anndataR)
library(reticulate)
reticulate::py_config()

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[which(args == "--input") + 1]
output_dir <- args[which(args == "--output") + 1]

if (length(input_dir) == 0 || length(output_dir) == 0) {
  stop("Missing --input or --output arguments")
}

# --- Construct dynamic paths ---
print("Loading fetal path")
fetal_path <- file.path(input_dir, "metaatlas_subsample.h5ad")
print(fetal_path)
print("Loading organoid path")
organoid_path <- file.path(input_dir, "rnh027_subsample.h5ad")
print(organoid_path)




print("Loading datasets...")
# Load the datasets directly as SingleCellExperiment objects
print("Loading fetal")
adata_fetal <- anndataR::read_h5ad(fetal_path)
print(adata_fetal)
print("Loading organoid")
adata_organoid <- anndataR::read_h5ad(organoid_path)
print(adata_organoid)


# Loading 


print("Converting fetal anndata to Single Cell Experiment")
sce_fetal <- SingleCellExperiment(
  assays = list(counts = t(adata_fetal$X)),
  colData = as.data.frame(adata_fetal$obs),
  rowData = as.data.frame(adata_fetal$var),
  # This is the exact correct way to pull obsm into a SingleCellExperiment!
  reducedDims = list(X_scVI = adata_fetal$obsm[["X_scVI"]]) 
)
print(sce_fetal)

print("Converting organoid anndata to Single Cell Experiment")
sce_organoid <- SingleCellExperiment(
  assays = list(counts = t(adata_organoid$X)),
  colData = as.data.frame(adata_organoid$obs),
  rowData = as.data.frame(adata_organoid$var),
  reducedDims = list(X_scVI = adata_organoid$obsm[["X_scVI"]])
)

print("Calculating logcounts for RIMA gene expression analysis...")
assays(sce_fetal)$logcounts <- log1p(assays(sce_fetal)$counts)
assays(sce_organoid)$logcounts <- log1p(assays(sce_organoid)$counts)

print(sce_organoid)
print(sce_fetal)


# --- NEW: Fix Dimensionality Reduction Names ---
print("Checking available dimensionality reductions...")
print(paste("Fetal available dims:", paste(reducedDimNames(sce_fetal), collapse=", ")))
print(paste("Organoid available dims:", paste(reducedDimNames(sce_organoid), collapse=", ")))


# Step 0: Define the neighbourhoods
print("Defining neighborhoods...")
define_neighbourhoods <- function(sce, prop_seeds, knn=10, reduced.dim="X_scVI"){
  n_components <- ncol(reducedDim(sce, reduced.dim))  
  mi <- Milo(sce)
  
  # Tell buildGraph to use X_scVI
  mi <- miloR::buildGraph(mi, k = knn, d = n_components, reduced.dim = reduced.dim)
  
  # Tell makeNhoods to ALSO use X_scVI (added reduced_dims argument here!)
  mi <- miloR::makeNhoods(mi, prop = prop_seeds, k = knn, d=n_components, reduced_dims = reduced.dim, refined = TRUE)
  
  return(mi)
}
mi_fetal <- define_neighbourhoods(sce_fetal, prop_seeds = 0.02)
mi_organoid <- define_neighbourhoods(sce_organoid, prop_seeds = 0.02)

# Step 1: Preprocess the Milo objects
print("Preprocessing Milo objects...")
milos <- preprocess_milos(mi_fetal, mi_organoid)

# --- DEBUG: Inspect Milo objects after preprocessing ---
for (i in seq_along(milos)) {
  print(paste("=== Post-preprocess Milo", i, "==="))
  print(paste("  class:", class(milos[[i]])))
  print(paste("  ncol (cells):", ncol(milos[[i]])))
  print(paste("  nrow (genes):", nrow(milos[[i]])))
  print(paste("  colData rows:", nrow(colData(milos[[i]]))))
  print(paste("  nhoods dim:", paste(dim(nhoods(milos[[i]])), collapse=" x ")))
  print(paste("  reducedDimNames:", paste(reducedDimNames(milos[[i]]), collapse=", ")))
  if ("X_scVI" %in% reducedDimNames(milos[[i]])) {
    print(paste("  X_scVI dim:", paste(dim(reducedDim(milos[[i]], "X_scVI")), collapse=" x ")))
  } else {
    print("  WARNING: X_scVI not found in reducedDims!")
    print(paste("  Available reducedDims:", paste(reducedDimNames(milos[[i]]), collapse=", ")))
  }
  # Check unified_celltype column
  ct <- colData(milos[[i]])$unified_celltype
  print(paste("  unified_celltype class:", class(ct)))
  print(paste("  unified_celltype length:", length(ct)))
  print(paste("  unified_celltype NAs:", sum(is.na(ct))))
  print(paste("  unified_celltype unique values:", paste(head(unique(as.character(ct)), 10), collapse=", ")))
  if (is.factor(ct)) {
    print(paste("  unified_celltype levels:", paste(head(levels(ct), 10), collapse=", ")))
  }
  # Check dataset column
  ds <- colData(milos[[i]])$dataset
  print(paste("  dataset class:", class(ds)))
  print(paste("  dataset unique values:", paste(head(unique(as.character(ds)), 10), collapse=", ")))
}

# Step 2: Calculate neighbourhood similarities
print("Calculating similarities...")
dt_sims <- calculate_similarities(milos, method = "spearman")

# --- NEW: Scrub the NAs out of your cell type column! (CORRECTED) ---
print("Patching missing cell types to prevent scrambling crash...")

# CRITICAL FIX: unified_celltype is a factor after preprocessing, so assigning
# a value like "Unknown" that isn't already a factor level silently produces NA
# instead of the value you wanted. We must convert to character first, patch
# the NAs/"nan"/empty strings, then store back as character.
for (i in seq_along(milos)) {
  ct <- as.character(colData(milos[[i]])$unified_celltype)
  n_na <- sum(is.na(ct) | ct == "nan" | ct == "")
  print(paste("  Milo", i, "- patching", n_na, "missing/nan/empty cell types"))
  ct[is.na(ct) | ct == "nan" | ct == ""] <- "Unknown"
  colData(milos[[i]])$unified_celltype <- ct
  print(paste("  Milo", i, "- unified_celltype NAs after patch:", sum(is.na(colData(milos[[i]])$unified_celltype))))
  print(paste("  Milo", i, "- unified_celltype class after patch:", class(colData(milos[[i]])$unified_celltype)))
}

# Step 3: Assess statistical significance of nhood-nhood similarity
print("Calculating significance...")
dt_sims_sig <- calculate_nhoodnhood_significance(
  milos, dt_sims,
  n_scrambles = 10,
  col_scramble_label = "unified_celltype",
  direction = "b"
)

# Step 4: Match significant nhood-nhood connections
print("Matching neighborhoods...")
print(paste("  Total rows in dt_sims_sig:", nrow(dt_sims_sig)))
print(paste("  Significant rows:", sum(dt_sims_sig$is_significant == TRUE)))
dt_match <- match_nhoods(dt_sims_sig[is_significant == TRUE])
print(paste("  Matched nhoods:", nrow(dt_match)))
print(paste("  dt_match columns:", paste(colnames(dt_match), collapse=", ")))

# Step 5: Downstream analysis (Top conserved genes)
print("Calculating conserved gene expression...")
dt_cope <- calculate_cope(milos, dt_match, genes = NULL)
dt_cope <- dt_cope[order(dt_cope$cope, na.last = FALSE), ]
print(paste("  COPE results:", nrow(dt_cope), "genes"))
print(paste("  Top 3 genes:", paste(tail(dt_cope$gene, 3), collapse=", ")))

# --- NEW: Save the actual results to your output directory! ---
print("Saving tables and plots to output directory...")

# 1. Save the data tables as CSVs
write.csv(dt_match, file.path(output_dir, "significant_matches.csv"), row.names=FALSE)
write.csv(dt_cope, file.path(output_dir, "top_conserved_genes.csv"), row.names=FALSE)

# 2. Save the plots as a multi-page PDF
# FIX: preprocess_milos can reorder/subset cells, which makes the nhood matrix
# dimensions inconsistent with colData. annotate_nhoods then crashes when it
# tries to set rownames. We wrap each plot in tryCatch so the script always
# finishes (CSVs are the critical output) and produces whatever plots succeed.
# --- DEBUG: Pre-plot diagnostics ---
print("=== Pre-plot diagnostics ===")
for (i in seq_along(milos)) {
  print(paste("Milo", i, "- ncol:", ncol(milos[[i]]),
              "nhoods:", ncol(nhoods(milos[[i]])),
              "nhood rows:", nrow(nhoods(milos[[i]]))))
  # Check if nhood matrix rows match number of cells
  if (nrow(nhoods(milos[[i]])) != ncol(milos[[i]])) {
    print(paste("  ** MISMATCH: nhood matrix has", nrow(nhoods(milos[[i]])),
                "rows but Milo has", ncol(milos[[i]]), "cells"))
  }
  # Check colnames consistency
  print(paste("  colData colnames:", paste(head(colnames(colData(milos[[i]])), 15), collapse=", ")))
  # Check if reducedDim rows match cells
  if ("X_scVI" %in% reducedDimNames(milos[[i]])) {
    rd_rows <- nrow(reducedDim(milos[[i]], "X_scVI"))
    print(paste("  X_scVI rows:", rd_rows, " cells:", ncol(milos[[i]])))
    if (rd_rows != ncol(milos[[i]])) {
      print("  ** MISMATCH: X_scVI rows != number of cells!")
    }
  }
}

pdf(file.path(output_dir, "RIMA_Analysis_Plots.pdf"), width = 12, height = 8)

tryCatch({
  print("Generating matches embedding plot...")
  print(plot_matches_embed(milos, dt_match,
        cols_color = c("unified_celltype", "unified_celltype"),
        dimred = "X_scVI"))
  print("  plot_matches_embed succeeded!")
}, error = function(e) {
  message("WARNING: plot_matches_embed failed: ", e$message)
  message("  Full traceback: ")
  message(paste(capture.output(traceback()), collapse="\n"))
})

tryCatch({
  print("Generating matches map plot...")
  print(plot_matches_map(milos, dt_match,
        cols_label = c("unified_celltype", "dataset")))
  print("  plot_matches_map succeeded!")
}, error = function(e) {
  message("WARNING: plot_matches_map failed: ", e$message)
  message("  Full traceback: ")
  message(paste(capture.output(traceback()), collapse="\n"))
})

tryCatch({
  print("Generating paired expression plot...")
  top_genes <- tail(dt_cope$gene, 3)
  print(paste("  Genes for paired expression:", paste(top_genes, collapse=", ")))
  if (length(top_genes) > 0) {
    print(plot_paired_expression(milos, dt_match, genes = top_genes))
    print("  plot_paired_expression succeeded!")
  } else {
    message("WARNING: No genes available for paired expression plot.")
  }
}, error = function(e) {
  message("WARNING: plot_paired_expression failed: ", e$message)
  message("  Full traceback: ")
  message(paste(capture.output(traceback()), collapse="\n"))
})

dev.off() # Close and save the PDF
# --- EXTRACTING INDIVIDUAL CELLS FROM MATCHED NEIGHBORHOODS ---
print("Extracting individual cells from matches...")

# Extract the cell-to-neighborhood sparse matrices from the Milo objects
fetal_nhoods <- nhoods(mi_fetal)
organoid_nhoods <- nhoods(mi_organoid)

# Create the extraction function
get_matched_cells <- function(match_row_index, matches_df, fetal_nhoods, organoid_nhoods) {
  
  # Get the neighborhood IDs from the first two columns of dt_match
  fetal_nhood_id <- as.character(matches_df[match_row_index, 1]) 
  organoid_nhood_id <- as.character(matches_df[match_row_index, 2])
  
  # Find which cells have a '1' for those specific neighborhoods
  fetal_cells <- rownames(fetal_nhoods)[fetal_nhoods[, fetal_nhood_id] == 1]
  organoid_cells <- rownames(organoid_nhoods)[organoid_nhoods[, organoid_nhood_id] == 1]
  
  return(list(
    fetal_matched_cells = fetal_cells,
    organoid_matched_cells = organoid_cells,
    similarity_score = matches_df$spearman[match_row_index] 
  ))
}

# Example: Get the cells for the very first matched pair
matched_pair_1 <- get_matched_cells(1, dt_match, fetal_nhoods, organoid_nhoods)

# --- NEW: INTERPRETATION AND SUMMARY ---
print("--- First Match Cell Summary ---")
print(paste("Similarity Score (Spearman):", matched_pair_1$similarity_score))
print(paste("Number of Fetal Cells:", length(matched_pair_1$fetal_matched_cells)))
print(paste("Number of Organoid Cells:", length(matched_pair_1$organoid_matched_cells)))

# Pull the metadata for these exact cells
fetal_meta <- colData(sce_fetal)[matched_pair_1$fetal_matched_cells, ]
organoid_meta <- colData(sce_organoid)[matched_pair_1$organoid_matched_cells, ]

print("--- Fetal Interpretation ---")
print("Cell Types:")
print(table(fetal_meta$unified_celltype))
print("Datasets:")
print(table(fetal_meta$dataset))

print("--- Organoid Interpretation ---")
print("Cell Types:")
print(table(organoid_meta$unified_celltype))
print("Datasets:")
print(table(organoid_meta$dataset))

# --- NEW: SAVE SUMMARIES TO A TEXT FILE ---
# This saves the tables to a text file so it uploads to your S3 bucket!
summary_file <- file.path(output_dir, "match_1_interpretation_summary.txt")
sink(summary_file)
cat("=== MATCH 1 INTERPRETATION SUMMARY ===\n")
cat("Spearman Similarity:", matched_pair_1$similarity_score, "\n\n")

cat("--- FETAL CELLS (n =", length(matched_pair_1$fetal_matched_cells), ") ---\n")
cat("Cell Types:\n")
print(table(fetal_meta$unified_celltype))
cat("\nDatasets:\n")
print(table(fetal_meta$dataset))

cat("\n--- ORGANOID CELLS (n =", length(matched_pair_1$organoid_matched_cells), ") ---\n")
cat("Cell Types:\n")
print(table(organoid_meta$unified_celltype))
cat("\nDatasets:\n")
print(table(organoid_meta$dataset))
sink()

# Save the raw barcodes to text files in the output directory
writeLines(matched_pair_1$fetal_matched_cells, file.path(output_dir, "match_1_fetal_cells.txt"))
writeLines(matched_pair_1$organoid_matched_cells, file.path(output_dir, "match_1_organoid_cells.txt"))

print("Pipeline complete!")