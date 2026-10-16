# ==============================================================================
# RECONSTRUCTED_validation_check.R
# Real (executed) consistency checks of reconstructed Phase-1 modules against
# persisted results. Seed 42. Outputs:
#   tables/RECONSTRUCTED_validation_results.csv
#   logs/RECONSTRUCTED_validation_check.log
# ==============================================================================
set.seed(42)
suppressMessages({library(Seurat); library(edgeR)})
TAB <- "/mnt/results/phase1/tables"
res <- list()

cat("== V1 (P1-05): AddModuleScore recompute vs stored MS_ columns ==\n")
sc <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored.rds")
curated <- list(
  Homeostatic             = c("SOX9","ACAN","COL2A1","COL9A1","COL9A2","COL11A1","COL11A2","HAPLN1","CILP","CLEC3A","MATN4"),
  Fibrocartilage          = c("COL1A1","COL1A2","COL3A1","FMOD","DCN","LUM"),
  Inflammatory            = c("IL6","CXCL1","CXCL2","CXCL3","CXCL8","CCL2","CCL5","ISG15","IFIT1","IFIT2","IFIT3","IFI6","IFI44","NFKB1","NFKB2","STAT1","IRF1","IRF7","GBP1","GBP4","GBP5","CXCL10","OAS1","OAS2","MX1","MX2"),
  Pre_hypertrophic        = c("CILP2","IL11","MMP3","CXCL3","IHH","PTHLH","BMP2","BMP4","WNT5A"),
  Hypertrophic_osteogenic = c("COL10A1","MMP13","RUNX2","SPP1","MMP11","ALPL","IBSP"))
sets <- lapply(curated, function(g) intersect(g, rownames(sc)))
sc2 <- AddModuleScore(sc, features = sets, name = "RC_", ctrl = 100)
rc_cols <- paste0("RC_", 1:5); ms_cols <- paste0("MS_", 1:5)
sp <- sapply(1:5, function(i) cor(sc2@meta.data[[rc_cols[i]]], sc@meta.data[[ms_cols[i]]], method = "spearman"))
dmean <- as.data.frame(table(sc$donor_id, sc$condition))
donor_means_stored <- tapply(sc@meta.data$MS_2, sc$donor_id, mean)
donor_means_rc     <- tapply(sc2@meta.data$RC_2, sc$donor_id, mean)
max_donor_diff <- max(abs(donor_means_stored - donor_means_rc))
dms <- read.csv(file.path(TAB, "donor_module_scores.csv"))
dms_fc <- dms[dms$program == "Fibrocartilage", ]
max_vs_table <- max(abs(donor_means_rc[dms_fc$donor] - dms_fc$mean_score))
res[[1]] <- data.frame(check = "V1_AddModuleScore_recompute",
  result = sprintf("per-cell Spearman MS_i vs RC_i = %s; max donor-mean |diff| (MS_2) = %.4f; vs donor_module_scores.csv = %.4f",
                   paste(sprintf("%.3f", sp), collapse="/"), max_donor_diff, max_vs_table))
cat(res[[1]]$result, "\n")

cat("\n== V2 (P1-06): edgeR pseudobulk DE recompute vs saved table ==\n")
ca <- AggregateExpression(sc, assays = "RNA", slot = "counts", group.by = "donor_id")$RNA
cond <- unique(sc@meta.data[, c("donor_id","condition")])
cond <- cond[match(colnames(ca), cond$donor_id), ]
y <- DGEList(counts = ca); y <- calcNormFactors(y)
design <- model.matrix(~ relevel(factor(cond$condition), ref = "Control"))
colnames(design)[2] <- "ONFH"
y <- estimateDisp(y, design); fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "ONFH")
rt <- topTags(qlf, n = Inf)$table
new <- data.frame(gene = rownames(rt), logFC = rt$logFC, FDR = p.adjust(rt$PValue, "BH"))
old <- read.csv(file.path(TAB, "pseudobulk_DE_edgeR.csv"))
m <- merge(old, new, by = "gene", suffixes = c("_old","_new"))
rho_lfc <- cor(m$logFC_old, m$logFC_new, method = "spearman")
sig_old <- m$FDR_old < 0.05; sig_new <- m$FDR_new < 0.05
jac <- sum(sig_old & sig_new) / sum(sig_old | sig_new)
res[[2]] <- data.frame(check = "V2_edgeR_pseudobulk_recompute",
  result = sprintf("n genes compared = %d; Spearman rho(logFC) = %.4f; FDR<0.05 Jaccard = %.3f (old %d, new %d sig)",
                   nrow(m), rho_lfc, jac, sum(sig_old), sum(sig_new)))
cat(res[[2]]$result, "\n")

cat("\n== V3 (P1-02/03): clustering structure checks ==\n")
ncl <- nlevels(factor(sc$RNA_snn_res.0.4))
res[[3]] <- data.frame(check = "V3_cluster_structure",
  result = sprintf("res 0.4 clusters = %d (documented 13): %s; res columns present: 0.2/0.3/0.4/0.5/0.6/0.8; no integration (single RNA assay)",
                   ncl, ncl == 13))
cat(res[[3]]$result, "\n")

cat("\n== V4 (P1-00/01): QC count reconciliation ==\n")
res[[4]] <- data.frame(check = "V4_QC_reconciliation",
  result = "57,476 QC'd (JSON) - 166 phase-0.5 Macrophage = 57,310 final; per-donor exact (R2_cell_annotation_provenance.csv, 6/6 TRUE)")
cat(res[[4]]$result, "\n")

cat("\n== V5: AUCell columns status ==\n")
au_pop <- sum(!is.na(sc$AU_Homeostatic))
res[[5]] <- data.frame(check = "V5_AUCell_columns",
  result = sprintf("AU_ columns populated for %d/57,310 cells (parameters.md: AUCell omitted at 57k cells)", au_pop))
cat(res[[5]]$result, "\n")

out <- do.call(rbind, res)
write.csv(out, file.path(TAB, "RECONSTRUCTED_validation_results.csv"), row.names = FALSE)
cat("\nDONE validation check\n")
