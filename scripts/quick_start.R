# Imports 
library(RIMA)
library(miloR)
library(SingleCellExperiment)
library(anndataR)
library(reticulate)
library(ggplot2) # <--- ADDED GGPLOT2 HERE
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

print("Converting fetal anndata to Single Cell Experiment")
sce_fetal <- SingleCellExperiment(
  assays = list(counts = t(adata_fetal$X)),
  colData = as.data.frame(adata_fetal$obs),
  rowData = as.data.frame(adata_fetal$var),
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
  
  mi <- miloR::buildGraph(mi, k = knn, d = n_components, reduced.dim = reduced.dim)
  mi <- miloR::makeNhoods(mi, prop = prop_seeds, k = knn, d=n_components, reduced_dims = reduced.dim, refined = TRUE)
  
  return(mi)
}
mi_fetal <- define_neighbourhoods(sce_fetal, prop_seeds = 0.005, knn = 30)
mi_organoid <- define_neighbourhoods(sce_organoid, prop_seeds = 0.005, knn = 30)

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
}

# Step 2: Calculate neighbourhood similarities
print("Calculating similarities...")
dt_sims <- calculate_similarities(milos, method = "spearman")

# --- NEW: Scrub the NAs out of your cell type column! (CORRECTED) ---
print("Patching missing cell types to prevent scrambling crash...")
for (i in seq_along(milos)) {
  ct <- as.character(colData(milos[[i]])$unified_celltype)
  n_na <- sum(is.na(ct) | ct == "nan" | ct == "")
  print(paste("  Milo", i, "- patching", n_na, "missing/nan/empty cell types"))
  ct[is.na(ct) | ct == "nan" | ct == ""] <- "Unknown"
  colData(milos[[i]])$unified_celltype <- ct
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
pdf(file.path(output_dir, "RIMA_Analysis_Plots.pdf"), width = 12, height = 8)

tryCatch({
  print("Generating matches embedding plot...")
  embed_plot <- plot_matches_embed(
    milos, dt_match,
    cols_color = c("unified_celltype", "unified_celltype"),
    dimred = "X_scVI"
  )
  embed_plot <- embed_plot + 
    labs(
      title = "RIMA Neighborhood Matching: Fetal Brain vs. SLUG-noid",
      subtitle = "Lines connect similar cellular states across the in vivo / in vitro barrier",
      x = "scVI Latent Space (Split by Dataset)",
      y = "scVI Latent Space"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray30")
    ) +
    annotate("text", x = -Inf, y = Inf, label = "Fetal Brain\n(In Vivo)", 
             hjust = -0.5, vjust = 1.5, size = 5, fontface = "bold") +
    annotate("text", x = Inf, y = Inf, label = "Brain Organoid\n(In Vitro)", 
             hjust = 1.5, vjust = 1.5, size = 5, fontface = "bold")

  print(embed_plot)
  print("  plot_matches_embed succeeded and customized!")
}, error = function(e) {
  message("WARNING: plot_matches_embed failed: ", e$message)
})

tryCatch({
  print("Generating matches map plot...")
  
  # 1. Save the plot to a variable
  match_heatmap <- plot_matches_map(milos, dt_match, cols_label = c("unified_celltype", "unified_celltype"))
  
  # 2. Overwrite the ugly default labels with clean ones
  match_heatmap <- match_heatmap + labs(
    title = "Neighborhood Matching Heatmap",
    x = "Fetal Brain Cell Types (In Vivo)",
    y = "SLUG-noid Cell Types (In Vitro)",
    fill = "Similarity Score" # Optional: cleans up the legend title too
  )
  
  # 3. Print it to the PDF
  print(match_heatmap)
  print("  plot_matches_map succeeded and customized!")
  
}, error = function(e) {
  message("WARNING: plot_matches_map failed: ", e$message)
})

tryCatch({
  print("Generating paired expression plot for SLUG-noid maturation genes...")
  
  # Define the specific genes the paper cares about (plus SCN2A as a positive control!)
  desired_genes <- c("MEF2C", "SATB2", "SCN2A")
  
  # Intersect with the actual genes in your dataset to prevent crashes if one got filtered out
  target_genes <- intersect(desired_genes, dt_cope$gene)
  
  if (length(target_genes) > 0) {
    # Generate the base plot
    paired_plot <- plot_paired_expression(milos, dt_match, genes = target_genes)
    
    # Add ggplot layers to clean up the axes and add a title
    paired_plot <- paired_plot + 
      labs(
        title = "Paired Gene Expression: Fetal vs Organoid",
        subtitle = "Assessing conservation of key neuronal maturation markers",
        x = "Fetal Expression (LogCounts)", 
        y = "Organoid Expression (LogCounts)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 12) # Makes the gene names bold!
      )
      
    print(paired_plot)
    print("  plot_paired_expression succeeded and customized!")
  } else {
    print("  WARNING: None of the target maturation genes were found in the dataset.")
  }
}, error = function(e) {
  message("WARNING: plot_paired_expression failed: ", e$message)
})

dev.off() 

# --- EXTRACTING INDIVIDUAL CELLS FROM MATCHED NEIGHBORHOODS ---
print("Extracting individual cells from matches...")
fetal_nhoods <- nhoods(mi_fetal)
organoid_nhoods <- nhoods(mi_organoid)

get_matched_cells <- function(match_row_index, matches_df, fetal_nhoods, organoid_nhoods) {
  fetal_nhood_id <- as.character(matches_df[match_row_index, 1]) 
  organoid_nhood_id <- as.character(matches_df[match_row_index, 2])
  fetal_cells <- rownames(fetal_nhoods)[fetal_nhoods[, fetal_nhood_id] == 1]
  organoid_cells <- rownames(organoid_nhoods)[organoid_nhoods[, organoid_nhood_id] == 1]
  
  return(list(
    fetal_matched_cells = fetal_cells,
    organoid_matched_cells = organoid_cells,
    similarity_score = matches_df$sim[match_row_index] 
  ))
}

matched_pair_1 <- get_matched_cells(1, dt_match, fetal_nhoods, organoid_nhoods)
fetal_meta <- colData(sce_fetal)[matched_pair_1$fetal_matched_cells, ]
organoid_meta <- colData(sce_organoid)[matched_pair_1$organoid_matched_cells, ]

summary_file <- file.path(output_dir, "match_1_interpretation_summary.txt")
sink(summary_file)
cat("=== MATCH 1 INTERPRETATION SUMMARY ===\n")
cat("Spearman Similarity:", matched_pair_1$similarity_score, "\n\n")
cat("--- FETAL CELLS (n =", length(matched_pair_1$fetal_matched_cells), ") ---\n")
print(table(fetal_meta$unified_celltype))
cat("\n--- ORGANOID CELLS (n =", length(matched_pair_1$organoid_matched_cells), ") ---\n")
print(table(organoid_meta$unified_celltype))
sink()

writeLines(matched_pair_1$fetal_matched_cells, file.path(output_dir, "match_1_fetal_cells.txt"))
writeLines(matched_pair_1$organoid_matched_cells, file.path(output_dir, "match_1_organoid_cells.txt"))

# ==========================================
# 1. CALCULATE AVERAGE FETAL AGE PER NEIGHBORHOOD (THE MISSING PIECE!)
# ==========================================
print("Calculating average gestational age for fetal neighborhoods...")
fetal_nhoods_mat <- nhoods(mi_fetal)
fetal_cell_ages <- colData(sce_fetal)$Gest_week_num

fetal_nhood_ages <- data.frame(
  fetal_nhood_id = colnames(fetal_nhoods_mat),
  avg_fetal_age = NA
)

for(i in 1:ncol(fetal_nhoods_mat)) {
  cells_in_nhood <- fetal_nhoods_mat[, i] == 1
  fetal_nhood_ages$avg_fetal_age[i] <- mean(fetal_cell_ages[cells_in_nhood], na.rm = TRUE)
}

# ==========================================
# 2. TRANSFER AGE AND ECM TO ORGANOIDS VIA RIMA
# ==========================================
print("Transferring ages, ECM conditions, and Cell Types to Organoid neighborhoods via RIMA...")
mapped_ages <- data.frame(
  fetal_nhood_id = as.character(dt_match[[1]]),
  organoid_nhood_id = as.character(dt_match[[2]]),
  similarity = dt_match$sim 
)

mapped_ages <- merge(mapped_ages, fetal_nhood_ages, by = "fetal_nhood_id", all.x = TRUE)

organoid_nhoods_mat <- nhoods(mi_organoid)
organoid_cell_days <- colData(sce_organoid)$day_num
organoid_cell_ecm <- as.character(colData(sce_organoid)$ecm) 
organoid_cell_types <- as.character(colData(sce_organoid)$unified_celltype) # <--- NEW

mapped_ages$organoid_culture_day <- NA
mapped_ages$ecm_condition <- NA 
mapped_ages$cell_type <- NA # <--- NEW

for(i in 1:nrow(mapped_ages)) {
  org_id <- mapped_ages$organoid_nhood_id[i]
  cells_in_nhood <- organoid_nhoods_mat[, org_id] == 1
  
  # Take the majority vote for the Day, completely eliminating "Fake Days"
  day_table <- table(organoid_cell_days[cells_in_nhood])
  mapped_ages$organoid_culture_day[i] <- as.numeric(names(day_table)[which.max(day_table)])
  
  # Take the majority vote for ECM
  ecm_table <- table(organoid_cell_ecm[cells_in_nhood])
  mapped_ages$ecm_condition[i] <- names(ecm_table)[which.max(ecm_table)]
  
  # Take the majority vote for the Cell Type
  type_table <- table(organoid_cell_types[cells_in_nhood]) # <--- NEW
  mapped_ages$cell_type[i] <- names(type_table)[which.max(type_table)] # <--- NEW
}

mapped_ages <- na.omit(mapped_ages)
mapped_ages$organoid_culture_day <- factor(mapped_ages$organoid_culture_day)

# ==========================================
# 3. PLOT THE RESULTS (SPLIT BY ECM)
# ==========================================
print("Generating Maturation Plot split by ECM and Cell Type...")
maturation_plot <- ggplot(mapped_ages, aes(x = organoid_culture_day, y = avg_fetal_age, fill = ecm_condition)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, color = "black", outlier.shape = NA, alpha = 0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
             aes(color = ecm_condition), size = 1.0, alpha = 0.6) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 3) + # <--- THIS IS THE MAGIC LINE
  labs(
    title = "Cell-Type Specific Maturation: SLUG-noids vs Controls",
    subtitle = "Tracking predicted in vivo age across developmental lineages",
    x = "Organoid Culture Day (In Vitro)",
    y = "Predicted Fetal Gestational Week (In Vivo)",
    fill = "ECM Condition",
    color = "ECM Condition"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 10), # Makes the cell type labels pop
    strip.background = element_rect(fill = "lightgray"),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "Organoid_Maturation_by_CellType.pdf"), plot = maturation_plot, width = 14, height = 10)
# ==========================================
# 4. TEMPORAL ALIGNMENT HEATMAP (DTW)
# ==========================================
print("Generating DTW Temporal Alignment Heatmap...")

# Ensure the dtw package is installed and loaded
if (!requireNamespace("dtw", quietly = TRUE)) {
  # FIX: Swapped to the secure HTTPS cloud mirror so the firewall doesn't block it
  install.packages("dtw", repos = "https://cloud.r-project.org")
}
library(dtw)

# 1. Map Time Metadata to ALL possible edges (not just the matched ones)
time_edges <- as.data.frame(dt_sims_sig)

# FIX: Force the first two columns to be named "id_1" and "id_2" so merge() doesn't crash
colnames(time_edges)[1:2] <- c("id_1", "id_2")

# Create a lookup table for Fetal Weeks (rounded to integers for the heatmap grid)
fetal_lookup <- data.frame(
  id_1 = fetal_nhood_ages$fetal_nhood_id,
  fetal_week = round(fetal_nhood_ages$avg_fetal_age)
)
# Create a lookup table for Fetal Weeks (rounded to integers for the heatmap grid)
fetal_lookup <- data.frame(
  id_1 = fetal_nhood_ages$fetal_nhood_id,
  fetal_week = round(fetal_nhood_ages$avg_fetal_age)
)

# Create a lookup table for Organoid Days
org_lookup <- data.frame(
  id_2 = colnames(organoid_nhoods_mat),
  org_day = NA
)
for(i in 1:nrow(org_lookup)) {
  cells_in_nhood <- organoid_nhoods_mat[, i] == 1
  day_table <- table(organoid_cell_days[cells_in_nhood])
  org_lookup$org_day[i] <- as.numeric(names(day_table)[which.max(day_table)])
}

# Merge lookups into the massive edge table
time_edges <- merge(time_edges, fetal_lookup, by="id_1")
time_edges <- merge(time_edges, org_lookup, by="id_2")
time_edges <- na.omit(time_edges)

# 2. Calculate the Fraction of Retained (Significant) Edges
edge_summary <- aggregate(is_significant ~ fetal_week + org_day, data = time_edges,
                          FUN = function(x) sum(x == TRUE) / length(x))
colnames(edge_summary)[3] <- "fraction_retained"

# 3. Create Distance Matrix for DTW
# DTW minimizes cost, so our cost is (1 - fraction_retained)
fet_weeks <- sort(unique(edge_summary$fetal_week))
org_days <- sort(unique(edge_summary$org_day))

cost_matrix <- matrix(1, nrow=length(fet_weeks), ncol=length(org_days),
                      dimnames=list(fet_weeks, org_days))

for(i in 1:nrow(edge_summary)) {
  r <- as.character(edge_summary$fetal_week[i])
  c <- as.character(edge_summary$org_day[i])
  cost_matrix[r, c] <- 1 - edge_summary$fraction_retained[i]
}

# 4. Calculate the DTW Path
# dtw() finds the optimal alignment through the cost matrix
dtw_res <- dtw(cost_matrix, keep=TRUE, step.pattern=symmetric2)

dtw_path <- data.frame(
  org_day = org_days[dtw_res$index2],
  fetal_week = fet_weeks[dtw_res$index1]
)

# 5. Plot with ggplot2
# Convert to factors so the grid plots cleanly
edge_summary$org_day <- factor(edge_summary$org_day, levels=org_days)
edge_summary$fetal_week <- factor(edge_summary$fetal_week, levels=fet_weeks)
dtw_path$org_day <- factor(dtw_path$org_day, levels=org_days)
dtw_path$fetal_week <- factor(dtw_path$fetal_week, levels=fet_weeks)

dtw_plot <- ggplot(edge_summary, aes(x=org_day, y=fetal_week)) +
  geom_tile(aes(fill=fraction_retained), color = "white") +
  scale_fill_viridis_c(option="cividis", name="Fraction\nRetained Edges") +
  # FIX: Added inherit.aes = FALSE so the arrow doesn't look for the fraction_retained column
  geom_path(data=dtw_path, aes(x=org_day, y=fetal_week, group=1), inherit.aes = FALSE,
            color="#FF6699", linewidth=2, arrow=arrow(length=unit(0.15, "inches"), type="closed")) +
  labs(
    title="DTW Maturation Alignment: SLUG-noids vs Fetal Brain",
    subtitle="The pink arrow tracks the optimal developmental path across...",
    x = "Organoid Age (Days)",   # <--- Added readable X axis label
    y = "Fetal Age (Weeks)"      # <--- Added readable Y axis label
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", size=16),
    axis.title = element_text(face="bold")
  )

ggsave(file.path(output_dir, "DTW_Alignment_Path.pdf"), plot = dtw_plot, width = 8, height = 7)
print("Pipeline complete! DTW plot saved successfully.")