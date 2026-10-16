#!/usr/bin/env R
# =============================================================================
# R2 Analysis 1: External cohort negative-control specificity analysis
# =============================================================================
# Purpose: Test the frozen Fibrocartilage signature in the two EXTERNAL bulk
# cohorts (GSE74089, GSE345158) against:
#   (1) 1000 random gene sets matched for gene number + mean expression
#   (2) Hallmark EMT, (3) generic ECM, (4) generic fibroblast signatures
# Method: ssGSEA (same scoring family as the original manuscript validation),
#         disease effect = mean(score_disease) - mean(score_control)
# Seeds: 42. No change to grouping, sample inclusion, or signature definition.
# =============================================================================

set.seed(42)
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages({
  library(GSVA); library(limma)
  library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)
  library(openxlsx)
})
N_PERM <- 1000
SEED <- 42
OUT_TAB <- "/mnt/results/phase1/tables"
OUT_FIG <- "/mnt/results/phase1/figures"
MSIGDB_DIR <- "/mnt/shared-workspace/shared/phase1/geo_data/msigdb"
MSIGDB_VER <- "MSigDB v2026.1 (Broad GMT)"

cat("=== R2 Analysis 1: External specificity ===\n")
cat("GSVA version:", as.character(packageVersion("GSVA")), "\n")
cat("Comparator source:", MSIGDB_VER, "\n\n")

# -----------------------------------------------------------------------------
# 1. Frozen Fibrocartilage signature (exact gene list from manuscript)
# -----------------------------------------------------------------------------
frozen <- read.csv(file.path(OUT_TAB, "frozen_signatures.csv"))
fibro_all <- frozen[frozen$signature_name == "Fibrocartilage", ]
cat("Frozen Fibrocartilage signature:", nrow(fibro_all), "genes\n")
print(fibro_all$gene)
# Expected direction per gene from the frozen table (all should be up_in_ONFH)
cat("Direction values:", paste(unique(fibro_all$direction), collapse = ","), "\n\n")

# -----------------------------------------------------------------------------
# 2. Build gene-level expression matrices (same preprocessing as manuscript)
# -----------------------------------------------------------------------------
# --- GSE74089: probe-level -> collapse highest-variance probe per symbol ---
expr74089_probe <- readRDS("/mnt/shared-workspace/shared/phase1/geo_data/GSE74089_expr.rds")
probe_annot <- readRDS("/mnt/shared-workspace/shared/phase1/geo_data/GSE74089_probe_annotation.rds")
annot <- probe_annot[probe_annot$probe %in% rownames(expr74089_probe) & !is.na(probe_annot$symbol), ]
expr_annot <- expr74089_probe[annot$probe, , drop = FALSE]
annot$variance <- apply(expr_annot, 1, var)
annot <- annot[order(annot$symbol, -annot$variance), ]
annot_dedup <- annot[!duplicated(annot$symbol), ]
X74089 <- as.matrix(expr74089_probe[annot_dedup$probe, , drop = FALSE])
rownames(X74089) <- annot_dedup$symbol
nfh_samples <- c("GSM1909503", "GSM1909504", "GSM1909506", "GSM1909507")
grp74089 <- ifelse(colnames(X74089) %in% nfh_samples, "disease", "control")
cat("GSE74089 gene-level matrix:", nrow(X74089), "genes x", ncol(X74089), "samples\n")

# --- GSE345158: TPM -> Ensembl->symbol -> collapse highest variance -> log2(TPM+1) ---
tpm <- read.xlsx("/mnt/shared-workspace/shared/phase1/geo_data/GSE345158_gene_tpm_exp.annotation.xlsx")
# xlsx stores all cells as text -> coerce sample columns to numeric
tpm_mat <- as.matrix(apply(tpm[, -1, drop = FALSE], 2, as.numeric))
rownames(tpm_mat) <- tpm$gene_id
colnames(tpm_mat) <- colnames(tpm)[-1]
cat("TPM non-numeric coercion NAs:", sum(is.na(tpm_mat)), "\n")
sym <- mapIds(org.Hs.eg.db, keys = rownames(tpm_mat), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
keep <- !is.na(sym)
tpm_mat <- tpm_mat[keep, ]
sym <- sym[keep]
# collapse duplicated symbols: highest variance on log2(TPM+1)
ltpm <- log2(tpm_mat + 1)
ord <- order(sym, -apply(ltpm, 1, var))
sym_ord <- sym[ord]; ltpm_ord <- ltpm[ord, ]
keep2 <- !duplicated(sym_ord)
X345158 <- ltpm_ord[keep2, ]
rownames(X345158) <- sym_ord[keep2]
sonfh_samples <- paste0("Case-", 1:5)
grp345158 <- ifelse(colnames(X345158) %in% sonfh_samples, "disease", "control")
cat("GSE345158 gene-level matrix:", nrow(X345158), "genes x", ncol(X345158), "samples\n")
de345_orig <- read.csv(file.path(OUT_TAB, "validation_GSE345158_DE.csv"))
cat("Original GSE345158 DE universe:", nrow(de345_orig),
    "genes; overlap with rebuilt matrix:", sum(rownames(X345158) %in% de345_orig$gene), "\n\n")

# -----------------------------------------------------------------------------
# 3. Cohort DE logFC (for directional concordance) - from manuscript tables
# -----------------------------------------------------------------------------
de74089 <- read.csv(file.path(OUT_TAB, "validation_GSE74089_DE.csv"))
lfc74089 <- de74089$logFC; names(lfc74089) <- de74089$gene
lfc345158 <- de345_orig$logFC; names(lfc345158) <- de345_orig$gene

# -----------------------------------------------------------------------------
# 4. Reproduction check: ssGSEA of all 8 frozen signatures vs original scores
# -----------------------------------------------------------------------------
sig_names <- unique(frozen$signature_name)
frozen_sets <- lapply(sig_names, function(s) frozen$gene[frozen$signature_name == s])
names(frozen_sets) <- sig_names

run_ssgsea <- function(X, sets, nthread = 4) {
  sets_f <- lapply(sets, function(g) unique(intersect(g, rownames(X))))
  sz <- sapply(sets_f, length)
  ok <- sz >= 2
  if (any(!ok)) cat("  [dropped sets with <2 genes present]:",
                    paste(names(sz)[!ok], collapse = ","), "\n")
  sets_f <- sets_f[ok]
  p <- ssgseaParam(X, sets_f, normalize = TRUE)
  sc <- gsva(p, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(nthread))
  t(sc)
}

cat("Reproduction check: ssGSEA of 8 frozen signatures...\n")
sc74089_frozen <- run_ssgsea(X74089, frozen_sets)
sc345158_frozen <- run_ssgsea(X345158, frozen_sets)

verify <- data.frame()
for (s in sig_names) {
  # GSE74089
  d7 <- mean(sc74089_frozen[grp74089 == "disease", s]) -
        mean(sc74089_frozen[grp74089 == "control", s])
  o7 <- read.csv(file.path(OUT_TAB, "validation_GSE74089_signature_scores.csv"))
  o7 <- o7$diff[o7$signature == s]
  # GSE345158
  d3 <- mean(sc345158_frozen[grp345158 == "disease", s]) -
        mean(sc345158_frozen[grp345158 == "control", s])
  o3 <- read.csv(file.path(OUT_TAB, "validation_GSE345158_signature_scores.csv"))
  o3 <- o3$diff[o3$signature == s]
  verify <- rbind(verify, data.frame(cohort = "GSE74089", signature = s,
      diff_rerun = d7, diff_original = o7, abs_diff = abs(d7 - o7)))
  verify <- rbind(verify, data.frame(cohort = "GSE345158", signature = s,
      diff_rerun = d3, diff_original = o3, abs_diff = abs(d3 - o3)))
}
cat("\n--- Reproduction check (ssGSEA disease effect, rerun vs original) ---\n")
print(verify, digits = 4)
write.csv(verify, file.path(OUT_TAB, "R2_extspec_verification.csv"), row.names = FALSE)

# Per-sample z-score check for Fibrocartilage (vs saved z-scores)
z7 <- as.vector(scale(sc74089_frozen[, "Fibrocartilage"])); names(z7) <- rownames(sc74089_frozen)
z3 <- as.vector(scale(sc345158_frozen[, "Fibrocartilage"])); names(z3) <- rownames(sc345158_frozen)
zs7 <- read.csv(file.path(OUT_TAB, "validation_GSE74089_zscores.csv"))
zs3 <- read.csv(file.path(OUT_TAB, "validation_GSE345158_zscores.csv"))
zs7 <- zs7[match(names(z7), zs7$sample), "Fibrocartilage"]
zs3 <- zs3[match(names(z3), zs3$sample), "Fibrocartilage"]
cat(sprintf("Z-score reproduction: GSE74089 rho=%.4f | GSE345158 rho=%.4f\n",
    cor(z7, zs7), cor(z3, zs3)))

# -----------------------------------------------------------------------------
# 5. Random gene sets: expression-bin matched (deciles), N_PERM sets, seed 42
# -----------------------------------------------------------------------------
make_random_sets <- function(X, sig_genes, n_perm = N_PERM, seed = SEED) {
  set.seed(seed)
  genes <- rownames(X)
  gm <- rowMeans(X)
  # decile bins of mean expression (rank-based: robust to tied values)
  bins <- cut(rank(gm), breaks = 10, labels = FALSE)
  names(bins) <- genes
  sig_in <- intersect(sig_genes, genes)
  sig_bins <- bins[sig_in]
  bin_counts <- table(sig_bins)
  sets <- vector("list", n_perm)
  for (i in seq_len(n_perm)) {
    g <- unlist(lapply(names(bin_counts), function(b) {
      cand <- setdiff(names(bins)[bins == as.integer(b)], names(sig_bins))
      if (length(cand) >= bin_counts[b]) sample(cand, bin_counts[b]) else
        sample(cand, bin_counts[b], replace = TRUE)
    }))
    sets[[i]] <- g
  }
  sets
}

# -----------------------------------------------------------------------------
# 6. Run ssGSEA: fibrocartilage + 1000 random + 3 comparators, per cohort
# -----------------------------------------------------------------------------
# Comparator gene lists from MSigDB v2026.1 official GMTs (Broad Institute).
# Note: the discovery-cohort comparison (R2 Task 4) used msigdbr 26.1.1; the
# GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION term is absent from MSigDB C5 GMT
# releases, so the generic ECM comparator here uses GOBP_EXTRACELLULAR_MATRIX_
# ASSEMBLY + GOBP_COLLAGEN_FIBRIL_ORGANIZATION (documented in Table 2).
cat("\nBuilding comparator gene lists from", MSIGDB_VER, "GMTs...\n")
parse_gmt_line <- function(line) strsplit(line, "\t", fixed = TRUE)[[1]][-(1:2)]
gmt_h <- readLines(file.path(MSIGDB_DIR, "h.all.v2026.1.Hs.symbols.gmt"))
emt_genes <- parse_gmt_line(grep("^HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION\t", gmt_h, value = TRUE))
gmt_c5 <- readLines(file.path(MSIGDB_DIR, "c5.go.bp.v2026.1.Hs.symbols.gmt"))
ecm_genes <- union(
  parse_gmt_line(grep("^GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY\t", gmt_c5, value = TRUE)),
  parse_gmt_line(grep("^GOBP_COLLAGEN_FIBRIL_ORGANIZATION\t", gmt_c5, value = TRUE)))
fib_lines <- grep("FIBROBLAST_PROLIFERATION", gmt_c5, value = TRUE)
fib_genes <- unique(unlist(lapply(fib_lines, parse_gmt_line)))
cat("  Hallmark EMT:", length(emt_genes), "| Generic ECM:", length(ecm_genes),
    "| Generic Fibroblast:", length(fib_genes), "\n")

comp_sets_full <- list(Hallmark_EMT = emt_genes, Generic_ECM = ecm_genes,
                       Generic_Fibroblast = fib_genes)

analyze_cohort <- function(X, grp, lfc, cohort_name, fibro_genes, seed = SEED) {
  cat("\n==== Cohort:", cohort_name, "====\n")
  sig_present <- intersect(fibro_genes, rownames(X))
  cat("Fibrocartilage genes present:", length(sig_present), "/", length(fibro_genes), "\n")

  rand_sets <- make_random_sets(X, sig_present, seed = seed)
  cat("Generated", length(rand_sets), "random sets (size",
      length(sig_present), ", expression-bin matched)\n")

  all_sets <- c(list(Fibrocartilage = sig_present), rand_sets, comp_sets_full)
  names(all_sets)[2:(N_PERM + 1)] <- paste0("RAND_", seq_len(N_PERM))
  t0 <- Sys.time()
  scores <- run_ssgsea(X, all_sets)
  cat("ssGSEA done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
      "min |", ncol(scores), "sets\n")

  is_d <- grp == "disease"; is_c <- grp == "control"
  diff_all <- colMeans(scores[is_d, , drop = FALSE]) - colMeans(scores[is_c, , drop = FALSE])

  # directional concordance: fraction of set genes whose cohort logFC sign matches
  # the FROZEN per-gene direction (23 up_in_ONFH + DCN up_in_Control)
  expected_up_fibro <- fibro_all$direction == "up_in_ONFH"
  names(expected_up_fibro) <- fibro_all$gene
  conc_one <- function(genes, expected_up = NULL) {
    g <- intersect(genes, names(lfc))
    if (length(g) == 0) return(NA)
    if (is.null(expected_up)) return(mean(lfc[g] > 0))  # random sets: null = P(up)
    mean((lfc[g] > 0) == expected_up[g])
  }
  conc_fibro <- conc_one(sig_present, expected_up = expected_up_fibro)

  obs_diff <- diff_all[["Fibrocartilage"]]
  null_diff <- diff_all[grep("^RAND_", names(diff_all))]
  # empirical p (two-sided, add-one)
  p_diff <- (sum(abs(null_diff) >= abs(obs_diff)) + 1) / (length(null_diff) + 1)
  percentile <- 100 * mean(null_diff < obs_diff)

  # concordance null: fraction of random-set genes up in disease (null branch:
  # no frozen direction -> P(logFC > 0), genome-wide baseline ~50%)
  conc_null <- sapply(rand_sets, conc_one)
  p_conc <- (sum(conc_null >= conc_fibro, na.rm = TRUE) + 1) / (sum(!is.na(conc_null)) + 1)

  # comparators
  comp_res <- data.frame()
  for (cn in names(comp_sets_full)) {
    g_present <- intersect(comp_sets_full[[cn]], rownames(X))
    cd <- diff_all[[cn]]
    tt <- t.test(scores[is_d, cn], scores[is_c, cn])
    se_diff <- sqrt(var(scores[is_d, cn]) / sum(is_d) + var(scores[is_c, cn]) / sum(is_c))
    comp_res <- rbind(comp_res, data.frame(
      cohort = cohort_name, signature = cn, n_genes = length(g_present),
      control_mean = mean(scores[is_c, cn]), disease_mean = mean(scores[is_d, cn]),
      disease_effect = cd, se_effect = se_diff, t_stat = unname(tt$statistic),
      p_value = tt$p.value,
      direction = ifelse(cd > 0, "up_in_disease", "down_in_disease"),
      concordance_rate = conc_one(g_present)))
  }

  # per-sample scores for fibrocartilage + comparators
  keep_cols <- c("Fibrocartilage", names(comp_sets_full))
  ps <- data.frame(sample = rownames(scores), group = grp, scores[, keep_cols],
                   check.names = FALSE)

  list(cohort = cohort_name, sig_present = sig_present,
       null_diff = null_diff, conc_null = conc_null,
       obs_diff = obs_diff, p_diff = p_diff, percentile = percentile,
       conc_fibro = conc_fibro, p_conc = p_conc, comp_res = comp_res,
       per_sample = ps, rand_sizes = sapply(rand_sets, length))
}

res74089 <- analyze_cohort(X74089, grp74089, lfc74089, "GSE74089", fibro_all$gene)
res345158 <- analyze_cohort(X345158, grp345158, lfc345158, "GSE345158", fibro_all$gene)

# -----------------------------------------------------------------------------
# 7. Tables
# -----------------------------------------------------------------------------
# Exact gene list used per cohort
genes_out <- do.call(rbind, lapply(list(res74089, res345158), function(r) {
  data.frame(cohort = r$cohort, gene = fibro_all$gene,
             in_cohort = fibro_all$gene %in% r$sig_present,
             selection_criterion = fibro_all$selection_criterion)
}))
write.csv(genes_out, file.path(OUT_TAB, "R2_extspec_fibrocartilage_genes.csv"), row.names = FALSE)

# Raw null distribution
raw <- do.call(rbind, lapply(list(res74089, res345158), function(r) {
  data.frame(cohort = r$cohort, set_id = seq_len(N_PERM),
             disease_effect = r$null_diff, concordance_rate = r$conc_null)
}))
write.csv(raw, file.path(OUT_TAB, "R2_extspec_random_null_raw.csv"), row.names = FALSE)

# Table 1: random signature analysis summary
tab1 <- do.call(rbind, lapply(list(res74089, res345158), function(r) {
  data.frame(cohort = r$cohort,
             n_random_sets = N_PERM, set_size = length(r$sig_present),
             observed_effect = r$obs_diff,
             null_mean = mean(r$null_diff), null_sd = sd(r$null_diff),
             empirical_p_effect = r$p_diff, percentile_effect = r$percentile,
             observed_concordance = r$conc_fibro,
             null_concordance_mean = mean(r$conc_null, na.rm = TRUE),
             empirical_p_concordance = r$p_conc)
}))
write.csv(tab1, file.path(OUT_TAB, "R2_extspec_random_results.csv"), row.names = FALSE)

# Table 2a: comparator gene lists (with source definition)
comp_def <- c(
  Hallmark_EMT = "MSigDB v2026.1 H: HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  Generic_ECM = "MSigDB v2026.1 C5 GO:BP: GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY union GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
  Generic_Fibroblast = "MSigDB v2026.1 C5 GO:BP: union of all GOBP_*FIBROBLAST_PROLIFERATION* sets")
comp_genes <- do.call(rbind, lapply(names(comp_sets_full), function(cn)
  data.frame(signature = cn, source = comp_def[[cn]],
             n_genes_full = length(comp_sets_full[[cn]]),
             gene = comp_sets_full[[cn]])))
write.csv(comp_genes, file.path(OUT_TAB, "R2_extspec_comparator_genes.csv"), row.names = FALSE)

# Table 2b: comparator ssGSEA results
comp_all <- rbind(res74089$comp_res, res345158$comp_res)
fibro_row <- do.call(rbind, lapply(list(res74089, res345158), function(r) {
  sc <- r$per_sample
  is_d <- sc$group == "disease"; is_c <- sc$group == "control"
  tt <- t.test(sc[is_d, "Fibrocartilage"], sc[is_c, "Fibrocartilage"])
  se_diff <- sqrt(var(sc[is_d, "Fibrocartilage"]) / sum(is_d) +
                    var(sc[is_c, "Fibrocartilage"]) / sum(is_c))
  data.frame(cohort = r$cohort, signature = "Fibrocartilage",
             n_genes = length(r$sig_present),
             control_mean = mean(sc[is_c, "Fibrocartilage"]),
             disease_mean = mean(sc[is_d, "Fibrocartilage"]),
             disease_effect = r$obs_diff, se_effect = se_diff,
             t_stat = unname(tt$statistic), p_value = tt$p.value,
             direction = "up_in_disease", concordance_rate = r$conc_fibro)
}))
tab2 <- rbind(fibro_row, comp_all)
write.csv(tab2, file.path(OUT_TAB, "R2_extspec_comparator_results.csv"), row.names = FALSE)

# Per-sample ssGSEA scores (fibrocartilage + comparators, both cohorts)
ps_all <- rbind(res74089$per_sample, res345158$per_sample)
write.csv(ps_all, file.path(OUT_TAB, "R2_extspec_ssgsea_scores.csv"), row.names = FALSE)

cat("\n--- Table 1: Random signature analysis ---\n"); print(tab1, digits = 3)
cat("\n--- Table 2: Comparator results ---\n"); print(tab2, digits = 3)

# -----------------------------------------------------------------------------
# 8. Figure: external specificity comparison (3 panels)
# -----------------------------------------------------------------------------
svg_ok <- requireNamespace("svglite", quietly = TRUE)
save_fig <- function(p, name, w = 10, h = 3.6) {
  if (svg_ok) svglite::svglite(file.path(OUT_FIG, paste0(name, ".svg")), width = w, height = h)
  else svg(file.path(OUT_FIG, paste0(name, ".svg")), width = w, height = h)
  print(p); dev.off()
  ggsave(file.path(OUT_FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

null_df <- raw
obs_df <- data.frame(cohort = c("GSE74089", "GSE345158"),
                     observed = c(res74089$obs_diff, res345158$obs_diff))

pA <- ggplot(null_df, aes(x = disease_effect, fill = cohort)) +
  geom_histogram(bins = 40, alpha = 0.55, color = "grey30", linewidth = 0.1) +
  geom_vline(data = obs_df, aes(xintercept = observed, color = cohort),
             linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~cohort, scales = "free_x") +
  scale_fill_manual(values = c(GSE74089 = "#9DBCD4", GSE345158 = "#E2C6B4"),
                    guide = "none") +
  scale_color_manual(values = c(GSE74089 = "#0279EE", GSE345158 = "#FF9400"),
                     guide = "none") +
  labs(x = "ssGSEA disease effect (random signatures)", y = "Random sets (n=1000)",
       title = "A") +
  theme_classic(base_family = "Liberation Sans") +
  theme(strip.background = element_blank())

conc_df <- raw
conc_obs <- data.frame(cohort = c("GSE74089", "GSE345158"),
                       observed = c(res74089$conc_fibro, res345158$conc_fibro))
pB <- ggplot(conc_df, aes(x = concordance_rate, fill = cohort)) +
  geom_histogram(bins = 40, alpha = 0.55, color = "grey30", linewidth = 0.1) +
  geom_vline(data = conc_obs, aes(xintercept = observed, color = cohort),
             linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~cohort, scales = "free_x") +
  scale_fill_manual(values = c(GSE74089 = "#9DBCD4", GSE345158 = "#E2C6B4"), guide = "none") +
  scale_color_manual(values = c(GSE74089 = "#0279EE", GSE345158 = "#FF9400"), guide = "none") +
  labs(x = "Directional concordance rate (random signatures)", y = "Random sets (n=1000)",
       title = "B") +
  theme_classic(base_family = "Liberation Sans") +
  theme(strip.background = element_blank())

tab2$signature <- factor(tab2$signature,
                         levels = c("Fibrocartilage", "Hallmark_EMT", "Generic_ECM", "Generic_Fibroblast"))
pC <- ggplot(tab2, aes(x = signature, y = disease_effect, fill = signature)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.62, color = "grey25", linewidth = 0.2) +
  geom_errorbar(aes(ymin = disease_effect - se_effect, ymax = disease_effect + se_effect),
                width = 0.18, linewidth = 0.3) +
  facet_wrap(~cohort, scales = "free_y") +
  scale_fill_manual(values = c(Fibrocartilage = "#0279EE", Hallmark_EMT = "#E9ED4C",
                               Generic_ECM = "#75A025", Generic_Fibroblast = "#FD9BED"),
                    guide = "none") +
  labs(x = NULL, y = "ssGSEA disease effect", title = "C") +
  theme_classic(base_family = "Liberation Sans") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(0, 8, 0, 0), size = 9))

library(patchwork)
fig <- (pA | pB) / pC + plot_layout(heights = c(1, 1.15))
save_fig(fig, "R2_Fig_extspec", w = 10.5, h = 7)

# -----------------------------------------------------------------------------
# 9. sessionInfo
# -----------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_TAB, "R2_extspec_sessionInfo.txt"))
cat("\n=== R2 Analysis 1 Complete ===\n")
