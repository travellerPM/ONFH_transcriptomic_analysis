###############################################################################
# 07_freeze_signatures.R
# Phase 1 Chunk 7: Freeze discovery signatures BEFORE examining validation data
# Part 8: Create frozen_chondrocyte_signatures.xlsx
###############################################################################
set.seed(42)
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(openxlsx)
})

cat("=== Freezing Discovery Signatures ===\n")
cat("NOTE: This is done BEFORE examining GSE345158 or GSE74089 validation data.\n\n")

# Load GSE311084 scored object for extracting program markers
obj <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored.rds")
if (is(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)

# Load edgeR results for effect sizes
edgeR_res <- read.csv("/mnt/results/phase1/tables/pseudobulk_DE_edgeR.csv")

# Program definitions (curated signatures)
programs <- list(
  Homeostatic = c("SOX9","ACAN","COL2A1","COL9A1","COL9A2","COL11A1","COL11A2",
                  "HAPLN1","CILP","CLEC3A","MATN4"),
  Fibrocartilage = c("COL1A1","COL1A2","COL3A1","FMOD","DCN","LUM"),
  Inflammatory = c("IL6","CXCL1","CXCL2","CXCL3","CXCL8","CCL2","CCL5",
                   "ISG15","IFIT1","IFIT2","IFIT3","IFI6","IFI44","NFKB1","NFKB2",
                   "STAT1","IRF1","IRF7","GBP1","GBP4","GBP5","CXCL10","OAS1","OAS2",
                   "MX1","MX2"),
  Pre_hypertrophic = c("CILP2","IL11","MMP3","CXCL3","IHH","PTHLH","BMP2","BMP4","WNT5A"),
  Hypertrophic_osteogenic = c("COL10A1","MMP13","RUNX2","SPP1","BGLAP","MMP11","ALPL","IBSP")
)

prog_names <- names(programs)

# For each program, extract:
# 1. Curated signature genes (with effect sizes from pseudobulk DE)
# 2. Top data-driven markers (from FindMarkers)
# 3. Direction (up in ONFH or up in Control)

# Get program markers (data-driven)
Idents(obj) <- "dominant_program"
program_markers <- list()
for (p in prog_names) {
  cat("  Extracting markers for", p, "...\n")
  mk <- tryCatch({
    FindMarkers(obj, ident.1 = p, test.use = "wilcox",
                min.pct = 0.10, logfc.threshold = 0.25, only.pos = TRUE,
                max.cells.per.ident = 200)
  }, error = function(e) NULL)
  if (!is.null(mk) && nrow(mk) > 0) {
    mk$gene <- rownames(mk)
    mk$program <- p
    program_markers[[p]] <- mk
  }
}

# Build frozen signature table
frozen <- data.frame()

for (p in prog_names) {
  # --- Curated genes ---
  curated <- programs[[p]]
  for (g in curated) {
    # Get effect size from edgeR if available
    er <- edgeR_res[edgeR_res$gene == g, ]
    lfc <- if (nrow(er) > 0) er$logFC else NA
    fdr <- if (nrow(er) > 0) er$FDR else NA
    direction <- if (!is.na(lfc)) ifelse(lfc > 0, "up_in_ONFH", "up_in_Control") else "NA"

    frozen <- rbind(frozen, data.frame(
      signature_name = p,
      gene = g,
      source_dataset = "GSE311084",
      state_program = p,
      direction = direction,
      selection_criterion = "curated_signature",
      logFC_pseudobulk = lfc,
      FDR_pseudobulk = fdr,
      stringsAsFactors = FALSE
    ))
  }

  # --- Top data-driven markers (top 20 by logFC, FDR < 0.05) ---
  if (!is.null(program_markers[[p]])) {
    mk <- program_markers[[p]]
    mk_sig <- mk[mk$p_val_adj < 0.05, ]
    mk_sig <- mk_sig[order(-mk_sig$avg_log2FC), ]
    top_mk <- head(mk_sig, 20)
    for (i in seq_len(nrow(top_mk))) {
      g <- top_mk$gene[i]
      if (!(g %in% curated)) {  # Don't duplicate curated genes
        er <- edgeR_res[edgeR_res$gene == g, ]
        lfc <- if (nrow(er) > 0) er$logFC else NA
        fdr <- if (nrow(er) > 0) er$FDR else NA
        direction <- if (!is.na(lfc)) ifelse(lfc > 0, "up_in_ONFH", "up_in_Control") else "NA"

        frozen <- rbind(frozen, data.frame(
          signature_name = p,
          gene = g,
          source_dataset = "GSE311084",
          state_program = p,
          direction = direction,
          selection_criterion = "data_driven_marker",
          logFC_pseudobulk = lfc,
          FDR_pseudobulk = fdr,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# --- Add composite signatures ---
# 6. Conserved ONFH remodeling program (programs present in both datasets)
# Based on cross-context analysis, all 5 programs are conserved
# The ONFH-associated programs (enriched in ONFH vs Control) are:
# Fibrocartilage (log2FC_prop=7.58), Homeostatic (4.65), Inflammatory (1.56)
conserved_onfh_genes <- unique(c(
  programs[["Fibrocartilage"]],
  programs[["Inflammatory"]],
  head(rownames(program_markers[["Fibrocartilage"]][program_markers[["Fibrocartilage"]]$p_val_adj < 0.05, ]), 15),
  head(rownames(program_markers[["Inflammatory"]][program_markers[["Inflammatory"]]$p_val_adj < 0.05, ]), 15)
))
for (g in conserved_onfh_genes) {
  er <- edgeR_res[edgeR_res$gene == g, ]
  lfc <- if (nrow(er) > 0) er$logFC else NA
  fdr <- if (nrow(er) > 0) er$FDR else NA
  frozen <- rbind(frozen, data.frame(
    signature_name = "Conserved_ONFH_remodeling",
    gene = g,
    source_dataset = "GSE311084",
    state_program = "Fibrocartilage+Inflammatory",
    direction = if (!is.na(lfc)) ifelse(lfc > 0, "up_in_ONFH", "up_in_Control") else "NA",
    selection_criterion = "conserved_disease_associated",
    logFC_pseudobulk = lfc,
    FDR_pseudobulk = fdr,
    stringsAsFactors = FALSE
  ))
}

# 7. Cartilage-context-specific program (GSE311084-enriched)
# Pre_hypertrophic is strongly Control-enriched in GSE311084
cartilage_specific_genes <- unique(c(
  programs[["Pre_hypertrophic"]],
  head(rownames(program_markers[["Pre_hypertrophic"]][program_markers[["Pre_hypertrophic"]]$p_val_adj < 0.05, ]), 15)
))
for (g in cartilage_specific_genes) {
  er <- edgeR_res[edgeR_res$gene == g, ]
  lfc <- if (nrow(er) > 0) er$logFC else NA
  fdr <- if (nrow(er) > 0) er$FDR else NA
  frozen <- rbind(frozen, data.frame(
    signature_name = "Cartilage_context_specific",
    gene = g,
    source_dataset = "GSE311084",
    state_program = "Pre_hypertrophic",
    direction = if (!is.na(lfc)) ifelse(lfc > 0, "up_in_ONFH", "up_in_Control") else "NA",
    selection_criterion = "cartilage_context_enriched",
    logFC_pseudobulk = lfc,
    FDR_pseudobulk = fdr,
    stringsAsFactors = FALSE
  ))
}

# 8. Cystic-lesion-specific program (GSE290411-enriched)
# Hypertrophic_osteogenic has lower concordance and is more represented in GSE290411
cystic_specific_genes <- unique(c(
  programs[["Hypertrophic_osteogenic"]],
  head(rownames(program_markers[["Hypertrophic_osteogenic"]][program_markers[["Hypertrophic_osteogenic"]]$p_val_adj < 0.05, ]), 15)
))
for (g in cystic_specific_genes) {
  er <- edgeR_res[edgeR_res$gene == g, ]
  lfc <- if (nrow(er) > 0) er$logFC else NA
  fdr <- if (nrow(er) > 0) er$FDR else NA
  frozen <- rbind(frozen, data.frame(
    signature_name = "Cystic_lesion_specific",
    gene = g,
    source_dataset = "GSE311084+GSE290411",
    state_program = "Hypertrophic_osteogenic",
    direction = if (!is.na(lfc)) ifelse(lfc > 0, "up_in_ONFH", "up_in_Control") else "NA",
    selection_criterion = "cystic_lesion_enriched",
    logFC_pseudobulk = lfc,
    FDR_pseudobulk = fdr,
    stringsAsFactors = FALSE
  ))
}

# Remove duplicates (keep first occurrence)
frozen <- frozen[!duplicated(frozen[, c("signature_name", "gene")]), ]

# Save
write.csv(frozen, "/mnt/results/phase1/tables/frozen_signatures.csv", row.names = FALSE)

# Create xlsx with multiple sheets
wb <- createWorkbook()
for (sig in unique(frozen$signature_name)) {
  sub <- frozen[frozen$signature_name == sig, ]
  sheet_name <- gsub("[^A-Za-z0-9_]", "_", substr(sig, 1, 30))
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, sub)
}
addWorksheet(wb, "All_Signatures")
writeData(wb, "All_Signatures", frozen)
saveWorkbook(wb, "/workspace/phase1/frozen_chondrocyte_signatures.xlsx", overwrite = TRUE)

cat("\n=== Frozen Signatures Summary ===\n")
cat("Total genes across all signatures:", nrow(frozen), "\n")
cat("Signatures:\n")
for (sig in unique(frozen$signature_name)) {
  cat("  ", sig, ":", sum(frozen$signature_name == sig), "genes\n")
}

rm(obj); gc()
cat("\n=== Signatures frozen. DO NOT alter based on validation results. ===\n")
