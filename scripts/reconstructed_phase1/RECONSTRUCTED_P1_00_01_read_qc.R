# ==============================================================================
# RECONSTRUCTED — NOT the original Phase-1 script.
# Original P1-00/P1-01 were lost in a workspace storage reset (README,
# "Provenance note"). Rebuilt 2026-10-15 from persisted artifacts:
#   - onfh_data/phase0_celltype_results.json  (per-sample raw/QC counts, medians)
#   - phase0_5/GSE311084_qc_summary.rds       (per-sample QC summary)
#   - EmptyDrops CR matrices under onfh_data/GSE311084/ and .../GSE290411/
# Parameters that could not be recovered are marked UNKNOWN.
# Validation: per-sample QC counts reconcile exactly with the JSON and with the
# final object (see R2_cell_annotation_provenance.csv); exact filter thresholds
# are UNKNOWN, so re-execution is UNVERIFIED.
# ==============================================================================
set.seed(42)
suppressMessages(library(Seurat))

BASE <- "/mnt/shared-workspace/shared/onfh_data"
samples_311084 <- c("GSM9316639_A-01", "GSM9316640_A-02", "GSM9316641_A-03",
                    "GSM9316642_A-06", "GSM9316643_A-07", "GSM9316644_A-08")

# ---- QC thresholds: UNKNOWN (not persisted). Placeholders below; the recorded
# ---- per-sample filtered counts are the authoritative reference.
QC_MIN_FEATURES <- 200   # UNKNOWN placeholder
QC_MIN_COUNTS   <- 500   # UNKNOWN placeholder
QC_MAX_MT       <- 15    # UNKNOWN placeholder (median pct.mt 2.4-4.9% post-QC)

read_qc_one <- function(dir_name) {
  m <- Read10x(file.path(BASE, "GSE311084", dir_name))
  o <- CreateSeuratObject(m, min.cells = 3, min.features = QC_MIN_FEATURES,
                          project = dir_name)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  o <- subset(o, subset = nCount_RNA > QC_MIN_COUNTS & percent.mt < QC_MAX_MT)
  o
}

objs <- lapply(samples_311084, read_qc_one)
names(objs) <- sub("GSM\\d+_", "", samples_311084)
merged <- merge(objs[[1]], objs[-1])
cat("GSE311084 reconstructed QC counts per sample:\n")
print(table(merged$orig.ident))

# Reference (persisted truth): n_cells_filtered per sample in
# phase0_celltype_results.json = 5892/9232/15838/10182/9458/6874 (total 57,476).
# GSE290411 follows the same pattern (9 samples; final cartilage-lineage subset
# = 2,113 cells; see RECONSTRUCTED_P1_04_annotation.R).
cat("\nNOTE: exact QC thresholds UNKNOWN -> full re-execution UNVERIFIED.\n")
