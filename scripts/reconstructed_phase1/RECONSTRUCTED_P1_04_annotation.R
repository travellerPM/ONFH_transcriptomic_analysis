# ==============================================================================
# RECONSTRUCTED — NOT the original Phase-1 script.
# Original P1-04 (phase-0.5 lineage re-annotation) lost. Rebuilt 2026-10-15.
# Evidence base (persisted):
#   - phase0_5/GSE311084_major_calls.rds : 11 integrated clusters -> cell_type
#     with rationale strings "cartilage_lineage (cart_n=5, fib_n=0, top_other=5)"
#     / "cartilage_fibrocartilage_candidate (cart_n=5, fib_n=5)" /
#     "Competitive: Macrophage has 4 supportive markers vs Chondrocyte 3".
#   - phase0_5/GSE311084_lineage.rds : per-cluster n_supportive_markers and
#     mean_marker_expression for 9 lineages (Chondrocyte, Fibroblast_stromal,
#     Endothelial, Monocyte, Macrophage, T_cell, Mast_cell, Osteoclast,
#     Neutrophil).
#   - phase0_celltype_results.json : 13-lineage marker panels used at phase 0
#     (fibroblast = COL1A1,COL1A2,COL3A1,DCN,LUM,PDGFRA; chondrocyte =
#     COL2A1,ACAN,SOX9,COL9A1,COL9A2,COL11A1,HAPLN1; etc.).
#   - Final per-cell lineage_call / cell_type in GSE311084_scored.rds.
# UNKNOWN: the exact phase-0.5 gene lists (only counts persist) and the
# supportive-marker detection threshold. The panels below are reconstructed
# from the phase-0 JSON. Full re-derivation UNVERIFIED; the per-cell
# lineage_call in the persisted object is the authoritative annotation.
# Provenance audit: R2_cell_identity_audit.md.
# ==============================================================================

panels_9lineage <- list(
  Chondrocyte        = c("COL2A1","ACAN","SOX9","COL9A1","COL9A2","COL11A1","HAPLN1"),
  Fibroblast_stromal = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"),
  Endothelial        = c("PECAM1","VWF","CDH5","CLDN5","EGFL7","FLT1"),
  Monocyte           = c("CD14","LYZ","FCN1","VCAN","FCGR3A"),
  Macrophage         = c("CD68","CD163","MRC1","C1QA","C1QB","C1QC","LST1","AIF1"),
  T_cell             = c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","IL7R","TRAC"),
  Mast_cell          = c("KIT","TPSAB1","CPA3","MS4A2"),
  Osteoclast         = c("CTSK","ACP5","CALCR","NFATC1"),
  Neutrophil         = c("S100A8","S100A9","S100A12","FCGR3B","CSF3R","FPR1")
)

# Rule (from rationale strings): per integrated cluster, count supportive
# markers per lineage (marker expressed above detection threshold in the
# cluster; exact threshold UNKNOWN); cartilage_lineage if Chondrocyte markers
# present and Fibroblast_stromal = 0; cartilage_fibrocartilage_candidate if
# both present; competitive non-cartilage call if another lineage dominates.
call_cluster <- function(supportive_counts) {
  cart <- supportive_counts[["Chondrocyte"]]
  fib  <- supportive_counts[["Fibroblast_stromal"]]
  others <- supportive_counts[setdiff(names(supportive_counts),
                                      c("Chondrocyte","Fibroblast_stromal"))]
  top_other <- max(others)
  if (top_other > cart & names(which.max(others)) != "Fibroblast_stromal") {
    return(paste0("Competitive: ", names(which.max(others)), " has ",
                  top_other, " supportive markers vs Chondrocyte ", cart))
  }
  if (cart > 0 & fib == 0) return("cartilage_lineage")
  if (cart > 0 & fib > 0)  return("cartilage_fibrocartilage_candidate")
  "unresolved"
}

cat("Reconstructed rule documented. Authoritative per-cell annotations:\n")
cat("  lineage_call + cell_type columns of GSE311084_scored.rds /\n")
cat("  GSE290411_scored.rds (all 2,113 cells cartilage_fibrocartilage_candidate).\n")
