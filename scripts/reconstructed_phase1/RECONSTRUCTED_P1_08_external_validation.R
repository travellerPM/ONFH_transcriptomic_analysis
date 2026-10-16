# ==============================================================================
# RECONSTRUCTED — NOT the original Phase-1 script.
# Original P1-08 lost. Rebuilt 2026-10-15.
# Evidence base:
#   - parameters.md: GSE74089 (Agilent 4x44K, GPL13497) probe->gene via Feature
#     Extraction + org.Hs.eg.db, 44,495 probes -> 20,401 genes, collapsing =
#     highest-variance probe per gene; limma DE; ssGSEA + mean z-score.
#     GSE345158: Ensembl->symbol, 18,834 genes, DESeq2; ssGSEA + mean z-score.
#   - Persisted inputs: geo_data/GSE74089_expr.rds, GSE74089_probe_annotation.rds,
#     GSE345158 series matrix + count/TPM xlsx.
#   - Persisted outputs: validation_GSE74089_zscores.csv,
#     validation_GSE345158_zscores.csv, validation_grading_summary.csv.
# Already re-validated in the revision phase by executable scripts:
#   - R2_task5_gse74089_rerun.R : independent DE re-run, Spearman rho(logFC)=1.0
#   - R2_ext_specificity.R reproduction check: ssGSEA disease effect rho=0.982
#     (GSE74089) / 0.787 (GSE345158); z-score reproduction same run.
# ==============================================================================
set.seed(42)
suppressMessages({library(GSVA); library(limma)})

GEO <- "/mnt/shared-workspace/shared/phase1/geo_data"
TAB <- "/mnt/results/phase1/tables"

run_ssgsea <- function(X, sets, nthread = 4) {
  sets_f <- lapply(sets, function(g) intersect(g, rownames(X)))
  sets_f <- sets_f[sapply(sets_f, length) >= 1]
  p <- ssgseaParam(X, sets_f, normalize = TRUE)
  gsva(p, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(nthread))
}

mean_z <- function(logfc, genes, sig) {          # per-signature mean z-score
  z <- scale(logfc); names(z) <- genes
  mean(z[intersect(sig, names(z))], na.rm = TRUE)
}

# ---- GSE74089 ----
X <- readRDS(file.path(GEO, "GSE74089_expr.rds"))
pa <- readRDS(file.path(GEO, "GSE74089_probe_annotation.rds"))
# probe -> symbol; collapse multi-mapping probes by highest variance
# design: NFH (4) vs Control (4); limma: lmFit/eBayes -> logFC per gene
# per-gene z-score of moderated logFC; per-signature mean z; ssGSEA per sample.
# ---- GSE345158 ----
# counts from GSE345158_gene_count_exp.annotation.xlsx; Ensembl->symbol via
# org.Hs.eg.db; DESeq2 (SONFH 5 vs Control 4); z-scores; ssGSEA.
cat("Reconstructed external-validation module defined.\n")
cat("DE and ssGSEA/z-score parts already re-validated (rho=1.0 / 0.982 / 0.787).\n")
