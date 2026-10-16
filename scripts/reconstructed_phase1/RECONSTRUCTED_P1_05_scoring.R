# ==============================================================================
# RECONSTRUCTED — NOT the original Phase-1 script.
# Original P1-05 lost. Rebuilt 2026-10-15.
# Evidence base: scored objects carry MS_1..MS_5 (AddModuleScore), *UC_
# (UCell) and AU_* (AUCell) columns for the 5 programs; parameters.md records
# "AddModuleScore + UCell (AUCell omitted: memory at 57k cells)" — the AU_
# columns in the persisted object are therefore expected to be unpopulated.
# Program gene sets = curated rows of tables/frozen_signatures.csv (frozen
# BEFORE validation; must not be altered).
# VALIDATED against saved results in RECONSTRUCTED_validation_check.R (V1, V5).
# ==============================================================================
set.seed(42)
suppressMessages({library(Seurat); library(UCell)})

TAB <- "/mnt/results/phase1/tables"
frozen <- read.csv(file.path(TAB, "frozen_signatures.csv"))

# Curated program gene sets (direction-agnostic for scoring; the 5 core programs)
curated <- list(
  Homeostatic            = c("SOX9","ACAN","COL2A1","COL9A1","COL9A2","COL11A1","COL11A2","HAPLN1","CILP","CLEC3A","MATN4"),
  Fibrocartilage         = c("COL1A1","COL1A2","COL3A1","FMOD","DCN","LUM"),
  Inflammatory           = c("IL6","CXCL1","CXCL2","CXCL3","CXCL8","CCL2","CCL5","ISG15","IFIT1","IFIT2","IFIT3","IFI6","IFI44","NFKB1","NFKB2","STAT1","IRF1","IRF7","GBP1","GBP4","GBP5","CXCL10","OAS1","OAS2","MX1","MX2"),
  Pre_hypertrophic       = c("CILP2","IL11","MMP3","CXCL3","IHH","PTHLH","BMP2","BMP4","WNT5A"),
  Hypertrophic_osteogenic= c("COL10A1","MMP13","RUNX2","SPP1","MMP11","ALPL","IBSP")
)

score_programs <- function(obj) {
  genes <- rownames(obj)
  sets <- lapply(curated, function(g) intersect(g, genes))
  # AddModuleScore (MS_1..MS_5 order as in the persisted object)
  obj <- AddModuleScore(obj, features = sets, name = "MS_", ctrl = 100)
  # UCell (rank-based; secondary scoring family)
  obj <- AddModuleScore_UCell(obj, features = sets, name = "_UC")
  obj
}

# NOTE: AddModuleScore samples control features stochastically; with seed 42
# donor-level means reproduce the persisted donor_module_scores.csv to within
# numerical noise (validated: see RECONSTRUCTED_validation_results.csv).
cat("Reconstructed scoring module defined (AddModuleScore + UCell).\n")
