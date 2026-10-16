# ==============================================================================
# RECONSTRUCTED — NOT the original Phase-1 script.
# Original P1-10 lost. Rebuilt 2026-10-15.
# Evidence base: persisted figures Fig1-Fig7 (SVG + PNG 300 dpi) under
# /mnt/results/phase1/figures/ and their source tables under tables/.
# Figure sources:
#   Fig1  UMAP overview            <- GSE311084_reclustered.rds (DimPlot)
#   Fig2  module scores            <- donor_module_scores.csv / scored object
#   Fig2b feature UMAPs            <- scored object
#   Fig3  volcano + abundance      <- pseudobulk_DE_edgeR.csv; composition
#   Fig4  pathway enrichment       <- GSEA/ORA tables
#   Fig5  TF activity              <- decoupleR/dorothea pseudobulk table
#   Fig6  validation summary       <- validation_grading_summary.csv
#   Fig7  robustness summary       <- robustness_*.csv
# All figures already exist and are preserved; this script documents their
# provenance for regeneration. ggplot2; svglite; PNG 300 dpi; Liberation Sans.
# ==============================================================================
suppressMessages(library(ggplot2))
cat("Figure provenance documented; persisted figures authoritative.\n")
