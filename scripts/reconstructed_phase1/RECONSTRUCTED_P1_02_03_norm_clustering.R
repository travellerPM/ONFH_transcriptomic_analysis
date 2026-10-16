# ==============================================================================
# RECONSTRUCTED — NOT the original Phase-1 script.
# Original P1-02/P1-03 lost in a workspace storage reset. Rebuilt 2026-10-15.
# Evidence base:
#   - reclustered objects contain a single RNA assay (log-normalized data +
#     scale.data) and pca/umap reductions -> standard LogNorm pipeline,
#     NO batch integration (also stated in resolution_selection_summary.md:
#     "expected with n=3 vs 3 without batch integration").
#   - Resolution sweep 0.2-0.8 (7 resolutions); GSE311084 res=0.4 selected
#     (13 clusters); GSE290411 res=0.3 (6 clusters).
#   - All res columns 0.2/0.3/0.4/0.5/0.6/0.8 persist in the object metadata.
# UNKNOWN: number of PCs, UMAP parameters, min.cells.pct for features.
# Full re-execution UNVERIFIED (hours of compute); structural checks validated
# in RECONSTRUCTED_validation_check.R (V3).
# ==============================================================================
set.seed(42)
suppressMessages(library(Seurat))

run_norm_clustering <- function(obj, res_selected, label) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)          # nPCs UNKNOWN
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  for (r in c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8)) {
    obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
    obj <- FindClusters(obj, resolution = r, verbose = FALSE)
  }
  cat(label, ": clusters at res", res_selected, "=",
      nlevels(factor(obj@meta.data[[paste0("RNA_snn_res.", res_selected)]])), "\n")
  obj
}

# rc311 <- run_norm_clustering(merged_qc_object, 0.4, "GSE311084")  # heavy; see note
# rc290 <- run_norm_clustering(gse290411_qc_object, 0.3, "GSE290411")
cat("Script documents the reconstructed pipeline; heavy re-execution UNVERIFIED.\n")
cat("Persisted objects GSE311084_reclustered.rds / GSE290411_reclustered.rds\n")
cat("remain the authoritative analysis inputs.\n")
