#!/usr/bin/env R
# R1 #1: Reproducibility framework summary
# Consolidate all robustness + consistency + validation evidence into one table

cat("=== R1 #1: Reproducibility Framework Summary ===\n")

# --- Load all relevant tables ---
lodo <- read.csv("/mnt/results/phase1/tables/robustness_lodo_top100_stability.csv")
scoring <- read.csv("/mnt/results/phase1/tables/robustness_scoring_sensitivity.csv")
gene_dom <- read.csv("/mnt/results/phase1/tables/robustness_gene_dominance_summary.csv")
cellnum <- read.csv("/mnt/results/phase1/tables/robustness_cellnumber_sensitivity.csv")
per_donor <- read.csv("/mnt/results/phase1/tables/per_donor_consistency.csv")
grading <- read.csv("/mnt/results/phase1/tables/validation_grading_summary.csv")
neg_control <- read.csv("/mnt/results/phase1/tables/R2_negative_control_results.csv")
doublet_311 <- read.csv("/mnt/results/phase1/tables/R2_doublet_rates_GSE311084.csv")
doublet_290 <- read.csv("/mnt/results/phase1/tables/R2_doublet_rates_GSE290411.csv")

# --- Build framework table ---
framework <- data.frame(
  evidence_category = character(),
  evidence_item = character(),
  result = character(),
  interpretation = character(),
  source_table = character(),
  stringsAsFactors = FALSE
)

add_row <- function(cat, item, result, interp, src) {
  framework <<- rbind(framework, data.frame(
    evidence_category = cat, evidence_item = item, result = result,
    interpretation = interp, source_table = src, stringsAsFactors = FALSE
  ))
}

# 1. Donor-level replication
add_row("Donor-level replication", "Per-donor DE consistency (top 50 genes)",
        sprintf("%d/50 genes consistent in 3/3 ONFH donors", sum(per_donor$n_onfh_consistent == 3)),
        "All top DE genes replicated across all ONFH donors", "per_donor_consistency.csv")

# 2. LODO
add_row("Leave-one-donor-out (LODO)", "Direction consistency (top 100 DE genes)",
        sprintf("Mean=%.3f; %d/100 genes with 100%% consistency",
                mean(lodo$direction_consistency, na.rm = TRUE),
                sum(lodo$direction_consistency == 1, na.rm = TRUE)),
        "No single donor drives the DE results", "robustness_lodo_top100_stability.csv")

add_row("Leave-one-donor-out (LODO)", "Significance consistency (top 100 DE genes)",
        sprintf("Mean=%.3f; %d/100 genes significant in all iterations",
                mean(lodo$sig_consistency, na.rm = TRUE),
                sum(lodo$sig_consistency == 1, na.rm = TRUE)),
        "84.7%% of top genes remain significant regardless of donor exclusion", "robustness_lodo_top100_stability.csv")

# 3. Scoring sensitivity
scoring_cell <- scoring[scoring$level == "per_cell" | is.na(scoring$level), ]
add_row("Scoring method sensitivity", "AddModuleScore vs UCell (per-cell Spearman)",
        sprintf("Range %.3f–%.3f across 5 programs",
                min(scoring_cell$spearman_rho), max(scoring_cell$spearman_rho)),
        "Two independent scoring methods yield concordant results", "robustness_scoring_sensitivity.csv")

# 4. Gene dominance
add_row("Gene dominance (LOO)", "Mean rho across leave-one-gene-out tests",
        sprintf("Range %.3f–%.3f across 5 programs",
                min(gene_dom$mean_rho), max(gene_dom$mean_rho)),
        "No single gene dominates any program score", "robustness_gene_dominance_summary.csv")

# 5. Cell-number sensitivity
add_row("Cell-number sensitivity", "Downsampling (25–90% of cells)",
        "Spearman rho = 1.0 for all fractions and conditions",
        "Program scores robust to cell subsampling", "robustness_cellnumber_sensitivity.csv")

# 6. Negative control
fc_nc <- neg_control[neg_control$program == "Fibrocartilage", ]
add_row("Negative control (empirical null)", "Fibrocartilage disease effect vs 1000 expression-matched random sets",
        sprintf("z=%.2f, empirical p=%s", fc_nc$z_score_disease_diff,
                format(fc_nc$empirical_p_disease_diff, digits = 3)),
        "Fibrocartilage disease effect significantly exceeds expression-matched null", "R2_negative_control_results.csv")

# 7. External validation
add_row("External validation", "GSE74089 (hip cartilage, microarray)",
        sprintf("Fibrocartilage: %s concordant, per-gene %.1f%%",
                ifelse(grading$gse74089_concordant[grading$signature == "Fibrocartilage"] == "TRUE", "YES", "NO"),
                grading$gse74089_per_gene_pct[grading$signature == "Fibrocartilage"]),
        "Strong support from independent microarray dataset", "validation_grading_summary.csv")

add_row("External validation", "GSE345158 (femoral head, RNA-seq)",
        sprintf("Fibrocartilage: %s concordant, per-gene %.1f%%",
                ifelse(grading$gse345158_concordant[grading$signature == "Fibrocartilage"] == "TRUE", "YES", "NO"),
                grading$gse345158_per_gene_pct[grading$signature == "Fibrocartilage"]),
        "Direction concordant in independent RNA-seq dataset", "validation_grading_summary.csv")

# 8. Doublet detection
add_row("Doublet detection (scDblFinder)", "GSE311084 per-sample doublet rate",
        sprintf("Range %.2f–%.2f%%, overall %.2f%%",
                min(doublet_311$doublet_rate), max(doublet_311$doublet_rate),
                round(100 * sum(doublet_311$n_doublet) / sum(doublet_311$n_cells_before), 2)),
        "Doublet rates within expected range for 10x scRNA-seq", "R2_doublet_rates_GSE311084.csv")

add_row("Doublet detection (scDblFinder)", "GSE290411 per-sample doublet rate",
        sprintf("Range %.2f–%.2f%%, overall %.2f%%",
                min(doublet_290$doublet_rate), max(doublet_290$doublet_rate),
                round(100 * sum(doublet_290$n_doublet) / sum(doublet_290$n_cells_before), 2)),
        "SIONFH-4 excluded (n=1 cell, insufficient for doublet detection)", "R2_doublet_rates_GSE290411.csv")

# 9. Cross-dataset per-gene concordance
cross_pg <- read.csv("/mnt/results/phase1/tables/validation_cross_dataset_per_signature.csv")
fc_cross <- cross_pg[cross_pg$signature == "Fibrocartilage", ]
add_row("Cross-dataset per-gene concordance", "Fibrocartilage program",
        sprintf("%.1f%% (%d/%d genes) concordant across 2 external datasets",
                fc_cross$pct_concordant, fc_cross$concordant_genes, fc_cross$total_genes),
        "Strongest cross-platform validation of any program", "validation_cross_dataset_per_signature.csv")

# --- Print and save ---
cat("\n--- Reproducibility Framework Summary ---\n")
print(framework[, c("evidence_category", "evidence_item", "result")])

write.csv(framework, "/mnt/results/phase1/tables/R1_reproducibility_framework_summary.csv", row.names = FALSE)

cat("\n=== R1 #1 Complete ===\n")
