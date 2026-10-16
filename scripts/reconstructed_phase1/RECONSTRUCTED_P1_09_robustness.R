# ==============================================================================
# RECONSTRUCTED — NOT the original Phase-1 script.
# Original P1-09 lost. Rebuilt 2026-10-15.
# Evidence base: persisted robustness tables (robustness_lodo_all_genes.csv,
# robustness_lodo_top100_stability.csv, robustness_cellnumber_sensitivity.csv,
# robustness_scoring_sensitivity.csv, robustness_gene_dominance*.csv,
# robustness_resolution_sensitivity.csv) and parameters.md:
#   - LODO: 6 iterations, edgeR pseudobulk each
#   - Cell-number sensitivity: downsampling 25-90% of cells
#   - Scoring sensitivity: AddModuleScore vs UCell concordance
#   - Gene dominance: leave-one-gene-out per program (10k-cell subset)
#   - Resolution sensitivity: program scores across 7 resolutions (0.2-0.8)
# Full re-execution is compute-heavy and UNVERIFIED; the persisted tables are
# the authoritative robustness results. This script documents the procedures
# with seed 42 for reproducibility of any future re-run.
# ==============================================================================
set.seed(42)

lodo_edger <- function(obj, donors) {
  # for each donor d: pseudobulk remaining 5 donors, edgeR ONFH vs Control,
  # record program-gene logFC/FDR; stability = rank correlation across folds
  NULL
}

downsample_scores <- function(obj, pcts = seq(25, 90, by = 5)) {
  # per pct: stratified cell downsample (seed 42), recompute AddModuleScore,
  # record donor-level program means; compare to full-object means
  NULL
}

leave_one_gene_out <- function(obj, n_subset = 10000) {
  # per program: drop each gene, recompute AddModuleScore on a stratified
  # 10k-cell subset, record donor-mean shift and program-call concordance
  NULL
}

resolution_sensitivity <- function(obj) {
  # program scores are resolution-independent (cell-level); the persisted
  # analysis records dominant_program composition across res 0.2-0.8
  NULL
}

cat("Reconstructed robustness procedures documented; persisted tables authoritative.\n")
