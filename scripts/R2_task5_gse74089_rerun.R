#!/usr/bin/env R
# R2 Task 5: GSE74089 probe collapsing re-run
# Re-run the probe-to-gene collapsing + limma DE from saved probe-level data
# Verify results are identical to the original run

suppressMessages({
  library(limma)
})

cat("=== R2 Task 5: GSE74089 Probe Collapsing Re-run ===\n")

# --- 1. Load saved probe-level data and annotation ---
expr <- readRDS("/mnt/shared-workspace/shared/phase1/geo_data/GSE74089_expr.rds")
probe_annot <- readRDS("/mnt/shared-workspace/shared/phase1/geo_data/GSE74089_probe_annotation.rds")

cat("Expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")
cat("Probe annotation:", nrow(probe_annot), "probes\n")

# --- 2. Re-run probe collapsing (identical method to original) ---
cat("\nRe-running probe collapsing...\n")

# Keep only probes with gene symbols
annot <- probe_annot[probe_annot$probe %in% rownames(expr) & !is.na(probe_annot$symbol), ]
cat("Probes with gene symbol:", nrow(annot), "\n")

# Collapse by highest variance (same as original)
expr_annot <- expr[annot$probe, , drop = FALSE]
probe_var <- apply(expr_annot, 1, var)
annot$variance <- probe_var[annot$probe]
annot <- annot[order(annot$symbol, -annot$variance), ]
annot_dedup <- annot[!duplicated(annot$symbol), ]
cat("Unique genes after dedup:", nrow(annot_dedup), "\n")

gene_expr <- expr[annot_dedup$probe, , drop = FALSE]
rownames(gene_expr) <- annot_dedup$symbol
cat("Gene-level matrix:", nrow(gene_expr), "x", ncol(gene_expr), "\n")

# --- 3. Re-run limma DE ---
cat("\nRe-running limma DE...\n")
nfh_samples <- c("GSM1909503", "GSM1909504", "GSM1909506", "GSM1909507")
group <- factor(ifelse(colnames(gene_expr) %in% nfh_samples, "NFH", "Control"),
                levels = c("Control", "NFH"))

design <- model.matrix(~group)
fit <- lmFit(gene_expr, design)
fit <- eBayes(fit)
de_rerun <- topTable(fit, coef = "groupNFH", number = Inf, sort.by = "P")
de_rerun$gene <- rownames(de_rerun)

cat("Genes tested:", nrow(de_rerun), "\n")
cat("Significant (adj.P < 0.05):", sum(de_rerun$adj.P.Val < 0.05), "\n")

# --- 4. Compare with original results ---
cat("\n--- Comparing with original results ---\n")
de_orig <- read.csv("/mnt/results/phase1/tables/validation_GSE74089_DE.csv")
cat("Original DE results:", nrow(de_orig), "genes\n")

# Merge and compare
merged <- merge(de_orig, de_rerun, by = "gene", suffixes = c("_orig", "_rerun"))
cat("Genes in both:", nrow(merged), "\n")

# Check logFC concordance
cor_logfc <- cor(merged$logFC_orig, merged$logFC_rerun, method = "spearman")
max_diff_logfc <- max(abs(merged$logFC_orig - merged$logFC_rerun), na.rm = TRUE)
cat(sprintf("Spearman rho (logFC orig vs rerun): %.6f\n", cor_logfc))
cat(sprintf("Max |logFC difference|: %.8f\n", max_diff_logfc))

# Check p-value concordance
cor_pval <- cor(merged$P.Value_orig, merged$P.Value_rerun, method = "spearman")
cat(sprintf("Spearman rho (P.Value orig vs rerun): %.6f\n", cor_pval))

# Check significant gene overlap
sig_orig <- de_orig$gene[de_orig$adj.P.Val < 0.05]
sig_rerun <- de_rerun$gene[de_rerun$adj.P.Val < 0.05]
overlap <- length(intersect(sig_orig, sig_rerun))
cat(sprintf("Significant genes: orig=%d, rerun=%d, overlap=%d (%.1f%%)\n",
    length(sig_orig), length(sig_rerun), overlap,
    100 * overlap / max(length(sig_orig), 1)))

# Check top 20 genes
top20_orig <- head(de_orig$gene, 20)
top20_rerun <- head(de_rerun$gene, 20)
cat("Top 20 genes identical:", all(top20_orig == top20_rerun), "\n")

# --- 5. Save verification results ---
verification <- data.frame(
  check = c("n_probes_input", "n_probes_with_symbol", "n_unique_genes",
             "n_genes_tested", "n_significant_adjP05",
             "spearman_rho_logFC", "max_abs_diff_logFC",
             "spearman_rho_pvalue", "n_sig_overlap", "top20_identical"),
  original = c(nrow(expr), nrow(annot), nrow(annot_dedup),
                nrow(de_orig), sum(de_orig$adj.P.Val < 0.05),
                NA, NA, NA, length(sig_orig), NA),
  rerun = c(nrow(expr), nrow(annot), nrow(annot_dedup),
             nrow(de_rerun), sum(de_rerun$adj.P.Val < 0.05),
             cor_logfc, max_diff_logfc, cor_pval, length(sig_rerun),
             all(top20_orig == top20_rerun)),
  stringsAsFactors = FALSE
)
cat("\n--- Verification Summary ---\n")
print(verification)

write.csv(verification, "/mnt/results/phase1/tables/R2_gse74089_rerun_verification.csv", row.names = FALSE)

# Save re-run DE results
write.csv(de_rerun, "/mnt/results/phase1/tables/R2_gse74089_DE_rerun.csv", row.names = FALSE)

cat("\n=== R2 Task 5 Complete ===\n")
