# ==============================================================================
# RECONSTRUCTED — NOT the original Phase-1 script.
# Original P1-06 lost. Rebuilt 2026-10-15.
# Evidence base: tables/pseudobulk_DE_edgeR.csv (18,450 genes; columns gene,
# logFC, logCPM, F, PValue, FDR), plus DESeq2 and limma-voom concordance
# tables; parameters.md: "edgeR primary; donor-level pseudobulk, n=3 vs 3".
# VALIDATED against saved results in RECONSTRUCTED_validation_check.R (V2).
# ==============================================================================
set.seed(42)
suppressMessages({library(Seurat); library(edgeR)})

RDS <- "/mnt/shared-workspace/shared/phase1"
TAB <- "/mnt/results/phase1/tables"

run_pseudobulk_edger <- function(obj, condition_col = "condition", donor_col = "donor_id") {
  # Sum counts per donor (pseudobulk; donor = biological unit)
  ca <- AggregateExpression(obj, assays = "RNA", layer = "counts",
                            group.by = donor_col)$RNA
  cond <- unique(obj@meta.data[, c(donor_col, condition_col)])
  cond <- cond[match(colnames(ca), cond[[donor_col]]), ]
  y <- DGEList(counts = ca)
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ relevel(factor(cond[[condition_col]]), ref = "Control"))
  colnames(design) <- c("(Intercept)", "ONFH")
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "ONFH")
  res <- topTags(qlf, n = Inf)$table
  data.frame(gene = rownames(res), logFC = res$logFC, logCPM = res$logCPM,
             F = res$F, PValue = res$PValue, FDR = p.adjust(res$PValue, "BH"))
}

# res <- run_pseudobulk_edger(readRDS(file.path(RDS, "GSE311084_reclustered.rds")))
# write.csv(res, file.path(TAB, "pseudobulk_DE_edgeR_RECONSTRUCTED_check.csv"), row.names = FALSE)
cat("Reconstructed edgeR pseudobulk module defined (TMM + QL F-test).\n")
cat("DESeq2 / limma-voom concordance variants follow the same pseudobulk input.\n")
