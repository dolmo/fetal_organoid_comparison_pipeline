# dtw_alignment.R
# Per-cell-type Dynamic Time Warping (DTW) alignment of the organoid culture-day
# trajectory onto fetal gestational weeks, using the fraction of retained
# significant edges as the alignment cost. One PDF is written per cell type.

run_dtw_alignment <- function(dt_sims_sig, fetal_nhood_ages, mi_organoid, sce_organoid, output_dir) {
  print("Generating DTW Temporal Alignment Heatmaps by Cell Type...")

  if (!requireNamespace("dtw", quietly = TRUE)) {
    install.packages("dtw", repos = "https://cloud.r-project.org")
  }
  library(dtw)

  organoid_nhoods_mat <- nhoods(mi_organoid)
  organoid_cell_days <- colData(sce_organoid)$day_num
  organoid_cell_types <- as.character(colData(sce_organoid)$unified_celltype)

  time_edges <- as.data.frame(dt_sims_sig)
  colnames(time_edges)[1:2] <- c("id_1", "id_2")

  fetal_lookup <- data.frame(
    id_1 = fetal_nhood_ages$fetal_nhood_id,
    fetal_week = round(fetal_nhood_ages$avg_fetal_age)
  )

  # Majority-vote culture day and cell type for each organoid neighbourhood.
  org_lookup <- data.frame(
    id_2 = colnames(organoid_nhoods_mat),
    org_day = NA,
    org_type = NA
  )
  for (i in 1:nrow(org_lookup)) {
    cells_in_nhood <- organoid_nhoods_mat[, i] == 1
    if (sum(cells_in_nhood) > 0) {
      day_table <- table(organoid_cell_days[cells_in_nhood])
      org_lookup$org_day[i] <- as.numeric(names(day_table)[which.max(day_table)])
      type_table <- table(organoid_cell_types[cells_in_nhood])
      org_lookup$org_type[i] <- names(type_table)[which.max(type_table)]
    }
  }

  time_edges <- merge(time_edges, fetal_lookup, by = "id_1")
  time_edges <- merge(time_edges, org_lookup, by = "id_2")
  time_edges <- na.omit(time_edges)

  unique_types <- unique(time_edges$org_type)

  for (target_type in unique_types) {
    print(paste("Running DTW for:", target_type))
    type_edges <- time_edges[time_edges$org_type == target_type, ]

    # Need enough timepoints in both dimensions to draw an alignment path;
    # skip lineages that only exist on a single day (e.g. Stressed cells).
    if (nrow(type_edges) < 5 || length(unique(type_edges$fetal_week)) < 2 ||
        length(unique(type_edges$org_day)) < 2) {
      print(paste("  Skipping", target_type, "- not enough timepoints for alignment."))
      next
    }

    # Fraction of retained significant edges per (fetal_week, org_day) cell.
    edge_summary <- aggregate(is_significant ~ fetal_week + org_day, data = type_edges,
                              FUN = function(x) sum(x == TRUE) / length(x))
    colnames(edge_summary)[3] <- "fraction_retained"

    fet_weeks <- sort(unique(edge_summary$fetal_week))
    org_days <- sort(unique(edge_summary$org_day))

    cost_matrix <- matrix(1, nrow = length(fet_weeks), ncol = length(org_days),
                          dimnames = list(fet_weeks, org_days))
    for (i in 1:nrow(edge_summary)) {
      r <- as.character(edge_summary$fetal_week[i])
      c <- as.character(edge_summary$org_day[i])
      cost_matrix[r, c] <- 1 - edge_summary$fraction_retained[i]
    }

    tryCatch({
      dtw_res <- dtw(cost_matrix, keep = TRUE, step.pattern = symmetric2)
      dtw_path <- data.frame(
        org_day = org_days[dtw_res$index2],
        fetal_week = fet_weeks[dtw_res$index1]
      )

      edge_summary$org_day <- factor(edge_summary$org_day, levels = org_days)
      edge_summary$fetal_week <- factor(edge_summary$fetal_week, levels = fet_weeks)
      dtw_path$org_day <- factor(dtw_path$org_day, levels = org_days)
      dtw_path$fetal_week <- factor(dtw_path$fetal_week, levels = fet_weeks)

      dtw_plot <- ggplot(edge_summary, aes(x = org_day, y = fetal_week)) +
        geom_tile(aes(fill = fraction_retained), color = "white") +
        scale_fill_viridis_c(option = "cividis", name = "Fraction\nRetained Edges") +
        geom_path(data = dtw_path, aes(x = org_day, y = fetal_week, group = 1), inherit.aes = FALSE,
                  color = "#FF6699", linewidth = 2,
                  arrow = arrow(length = unit(0.15, "inches"), type = "closed")) +
        labs(
          title = paste("DTW Maturation:", target_type),
          subtitle = "Tracking the optimal developmental path for this specific lineage",
          x = "Organoid Age (Days)",
          y = "Fetal Age (Weeks)"
        ) +
        theme_minimal() +
        theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 16))

      safe_name <- gsub(" ", "_", target_type)
      file_name <- paste0("DTW_Alignment_", safe_name, ".pdf")
      ggsave(file.path(output_dir, file_name), plot = dtw_plot, width = 8, height = 7)
      print(paste("  Saved:", file_name))
    }, error = function(e) {
      print(paste("  DTW calculation failed for", target_type, ":", e$message))
    })
  }
}
