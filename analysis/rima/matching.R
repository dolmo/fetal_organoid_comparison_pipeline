# matching.R
# Neighbourhood-to-neighbourhood similarity significance, matching, and
# conserved-gene (COPE) scoring using the RIMA framework.

# Assess statistical significance of nhood-nhood similarity via label scrambling.
calculate_significance <- function(milos, dt_sims, n_scrambles = 10,
                                   col_scramble_label = "unified_celltype",
                                   direction = "b") {
  print("Calculating significance...")
  calculate_nhoodnhood_significance(
    milos, dt_sims,
    n_scrambles = n_scrambles,
    col_scramble_label = col_scramble_label,
    direction = direction
  )
}

# Match the significant nhood-nhood connections.
match_significant_nhoods <- function(dt_sims_sig) {
  print("Matching neighborhoods...")
  print(paste("  Total rows in dt_sims_sig:", nrow(dt_sims_sig)))
  print(paste("  Significant rows:", sum(dt_sims_sig$is_significant == TRUE)))
  dt_match <- match_nhoods(dt_sims_sig[is_significant == TRUE])
  print(paste("  Matched nhoods:", nrow(dt_match)))
  dt_match
}

# Calculate conserved gene expression (COPE) across matched neighbourhoods,
# sorted ascending by score.
calculate_conserved_genes <- function(milos, dt_match, genes = NULL) {
  print("Calculating conserved gene expression...")
  dt_cope <- calculate_cope(milos, dt_match, genes = genes)
  dt_cope <- dt_cope[order(dt_cope$cope, na.last = FALSE), ]
  print(paste("  COPE results:", nrow(dt_cope), "genes"))
  print(paste("  Top 3 genes:", paste(tail(dt_cope$gene, 3), collapse = ", ")))
  dt_cope
}
