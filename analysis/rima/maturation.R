# maturation.R
# Transfer the average fetal gestational age onto organoid neighbourhoods via
# the RIMA matches, then plot cell-type-specific maturation by ECM condition.

# Average gestational week of the cells in each fetal neighbourhood.
compute_fetal_nhood_ages <- function(mi_fetal, sce_fetal) {
  print("Calculating average gestational age for fetal neighborhoods...")
  fetal_nhoods_mat <- nhoods(mi_fetal)
  fetal_cell_ages <- colData(sce_fetal)$Gest_week_num

  fetal_nhood_ages <- data.frame(
    fetal_nhood_id = colnames(fetal_nhoods_mat),
    avg_fetal_age = NA
  )
  for (i in 1:ncol(fetal_nhoods_mat)) {
    cells_in_nhood <- fetal_nhoods_mat[, i] == 1
    fetal_nhood_ages$avg_fetal_age[i] <- mean(fetal_cell_ages[cells_in_nhood], na.rm = TRUE)
  }
  fetal_nhood_ages
}

# Map fetal age onto each matched organoid neighbourhood, plus a majority-vote
# culture day, ECM condition, and cell type per organoid neighbourhood.
build_mapped_ages <- function(dt_match, fetal_nhood_ages, mi_organoid, sce_organoid) {
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
  organoid_cell_types <- as.character(colData(sce_organoid)$unified_celltype)

  mapped_ages$organoid_culture_day <- NA
  mapped_ages$ecm_condition <- NA
  mapped_ages$cell_type <- NA

  for (i in 1:nrow(mapped_ages)) {
    org_id <- mapped_ages$organoid_nhood_id[i]
    cells_in_nhood <- organoid_nhoods_mat[, org_id] == 1

    # Majority vote for the culture day (eliminates spurious "fake days").
    day_table <- table(organoid_cell_days[cells_in_nhood])
    mapped_ages$organoid_culture_day[i] <- as.numeric(names(day_table)[which.max(day_table)])

    # Majority vote for the ECM condition.
    ecm_table <- table(organoid_cell_ecm[cells_in_nhood])
    mapped_ages$ecm_condition[i] <- names(ecm_table)[which.max(ecm_table)]

    # Majority vote for the cell type.
    type_table <- table(organoid_cell_types[cells_in_nhood])
    mapped_ages$cell_type[i] <- names(type_table)[which.max(type_table)]
  }

  mapped_ages <- na.omit(mapped_ages)
  mapped_ages$organoid_culture_day <- factor(mapped_ages$organoid_culture_day)
  mapped_ages
}

# Boxplot of predicted fetal age vs organoid culture day, split by ECM condition
# and faceted by cell type.
plot_maturation_by_celltype <- function(mapped_ages) {
  print("Generating Maturation Plot split by ECM and Cell Type...")
  ggplot(mapped_ages, aes(x = organoid_culture_day, y = avg_fetal_age, fill = ecm_condition)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, color = "black",
                 outlier.shape = NA, alpha = 0.8) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
               aes(color = ecm_condition), size = 1.0, alpha = 0.6) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
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
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "lightgray"),
      legend.position = "bottom"
    )
}
