# cell_extraction.R
# Pull the individual cells out of a matched neighbourhood pair and write a
# human-readable interpretation summary.

# Return the fetal and organoid cell barcodes belonging to a matched pair.
get_matched_cells <- function(match_row_index, matches_df, fetal_nhoods, organoid_nhoods) {
  fetal_nhood_id <- as.character(matches_df[match_row_index, 1])
  organoid_nhood_id <- as.character(matches_df[match_row_index, 2])
  fetal_cells <- rownames(fetal_nhoods)[fetal_nhoods[, fetal_nhood_id] == 1]
  organoid_cells <- rownames(organoid_nhoods)[organoid_nhoods[, organoid_nhood_id] == 1]

  list(
    fetal_matched_cells = fetal_cells,
    organoid_matched_cells = organoid_cells,
    similarity_score = matches_df$sim[match_row_index]
  )
}

# Write the cell barcodes and a cell-type breakdown for one matched pair.
write_match_summary <- function(match_index, dt_match, mi_fetal, mi_organoid,
                                sce_fetal, sce_organoid, output_dir) {
  print("Extracting individual cells from matches...")
  fetal_nhoods <- nhoods(mi_fetal)
  organoid_nhoods <- nhoods(mi_organoid)

  matched_pair <- get_matched_cells(match_index, dt_match, fetal_nhoods, organoid_nhoods)
  fetal_meta <- colData(sce_fetal)[matched_pair$fetal_matched_cells, ]
  organoid_meta <- colData(sce_organoid)[matched_pair$organoid_matched_cells, ]

  summary_file <- file.path(output_dir, paste0("match_", match_index, "_interpretation_summary.txt"))
  sink(summary_file)
  cat("=== MATCH", match_index, "INTERPRETATION SUMMARY ===\n")
  cat("Spearman Similarity:", matched_pair$similarity_score, "\n\n")
  cat("--- FETAL CELLS (n =", length(matched_pair$fetal_matched_cells), ") ---\n")
  print(table(fetal_meta$unified_celltype))
  cat("\n--- ORGANOID CELLS (n =", length(matched_pair$organoid_matched_cells), ") ---\n")
  print(table(organoid_meta$unified_celltype))
  sink()

  writeLines(matched_pair$fetal_matched_cells,
             file.path(output_dir, paste0("match_", match_index, "_fetal_cells.txt")))
  writeLines(matched_pair$organoid_matched_cells,
             file.path(output_dir, paste0("match_", match_index, "_organoid_cells.txt")))
}
