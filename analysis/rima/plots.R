# plots.R
# RIMA analysis figures, collected into a single multi-page PDF.

# Embedding plot connecting matched fetal / organoid neighbourhoods.
plot_matches_embedding <- function(milos, dt_match) {
  embed_plot <- plot_matches_embed(
    milos, dt_match,
    cols_color = c("unified_celltype", "unified_celltype"),
    dimred = "X_scVI"
  )
  embed_plot +
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
}

# Heatmap of the mean Spearman similarity between matched fetal and organoid
# cell types (majority cell type per neighbourhood).
plot_similarity_heatmap <- function(dt_match, sce_fetal, sce_organoid, mi_fetal, mi_organoid) {
  fetal_nhoods_mat <- nhoods(mi_fetal)
  organoid_nhoods_mat <- nhoods(mi_organoid)

  # Assign each fetal neighbourhood its majority cell type.
  fetal_cell_types <- as.character(colData(sce_fetal)$unified_celltype)
  fetal_nhood_types <- data.frame(fetal_nhood_id = colnames(fetal_nhoods_mat), fetal_type = NA)
  for (i in 1:ncol(fetal_nhoods_mat)) {
    cells_in_nhood <- fetal_nhoods_mat[, i] == 1
    t_tab <- table(fetal_cell_types[cells_in_nhood])
    fetal_nhood_types$fetal_type[i] <- names(t_tab)[which.max(t_tab)]
  }

  # Assign each organoid neighbourhood its majority cell type.
  organoid_cell_types <- as.character(colData(sce_organoid)$unified_celltype)
  org_nhood_types <- data.frame(org_nhood_id = colnames(organoid_nhoods_mat), org_type = NA)
  for (i in 1:ncol(organoid_nhoods_mat)) {
    cells_in_nhood <- organoid_nhoods_mat[, i] == 1
    t_tab <- table(organoid_cell_types[cells_in_nhood])
    org_nhood_types$org_type[i] <- names(t_tab)[which.max(t_tab)]
  }

  # Merge similarity scores with their assigned cell types and average per pair.
  sim_df <- data.frame(
    fetal_nhood_id = as.character(dt_match[[1]]),
    org_nhood_id = as.character(dt_match[[2]]),
    sim = dt_match$sim
  )
  sim_df <- merge(sim_df, fetal_nhood_types, by = "fetal_nhood_id")
  sim_df <- merge(sim_df, org_nhood_types, by.x = "org_nhood_id", by.y = "org_nhood_id")
  sim_agg <- aggregate(sim ~ fetal_type + org_type, data = sim_df, FUN = mean)

  ggplot(sim_agg, aes(x = fetal_type, y = org_type, fill = sim)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "plasma", name = "Mean\nSimilarity") +
    labs(
      title = "True Transcriptomic Similarity: Fetal vs Organoid",
      subtitle = "Color intensity represents the average Spearman similarity between mapped cell states",
      x = "Fetal Brain Cell Types (In Vivo)",
      y = "SLUG-noid Cell Types (In Vitro)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )
}

# Paired fetal-vs-organoid expression of key neuronal maturation markers.
# Returns NULL if none of the target genes are present.
plot_paired_maturation_genes <- function(milos, dt_match, dt_cope,
                                         desired_genes = c("MEF2C", "SATB2", "SCN2A")) {
  # Intersect with genes actually present so a filtered-out gene cannot crash us.
  target_genes <- intersect(desired_genes, dt_cope$gene)
  if (length(target_genes) == 0) {
    print("  WARNING: None of the target maturation genes were found in the dataset.")
    return(NULL)
  }

  paired_plot <- plot_paired_expression(milos, dt_match, genes = target_genes)
  paired_plot +
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
      strip.text = element_text(face = "bold", size = 12)
    )
}

# Render all three figures into one multi-page PDF. Each plot is wrapped in
# tryCatch so a single failure does not abort the others.
save_analysis_plots <- function(milos, dt_match, dt_cope,
                                sce_fetal, sce_organoid, mi_fetal, mi_organoid,
                                pdf_path) {
  print("Saving plots to output directory...")
  pdf(pdf_path, width = 12, height = 8)
  on.exit(dev.off(), add = TRUE)

  tryCatch({
    print("Generating matches embedding plot...")
    print(plot_matches_embedding(milos, dt_match))
  }, error = function(e) message("WARNING: plot_matches_embed failed: ", e$message))

  tryCatch({
    print("Generating Custom Average Similarity Heatmap...")
    print(plot_similarity_heatmap(dt_match, sce_fetal, sce_organoid, mi_fetal, mi_organoid))
  }, error = function(e) message("WARNING: Custom Similarity Heatmap failed: ", e$message))

  tryCatch({
    print("Generating paired expression plot for SLUG-noid maturation genes...")
    paired_plot <- plot_paired_maturation_genes(milos, dt_match, dt_cope)
    if (!is.null(paired_plot)) print(paired_plot)
  }, error = function(e) message("WARNING: plot_paired_expression failed: ", e$message))
}
