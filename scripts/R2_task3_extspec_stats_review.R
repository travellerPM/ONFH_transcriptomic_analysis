# ==============================================================================
# R2_task3_extspec_stats_review.R
# Reviewer-2 integrity audit, Task 3: statistical review of the external-cohort
# negative-control analysis (no re-running of ssGSEA; no change to results in
# place — original and corrected statistics are reported side by side).
#
# 3.1 Directional concordance
#   Observed statistic  : fraction of the 23 present Fibrocartilage genes whose
#                         cohort logFC sign matches the FROZEN per-gene direction
#                         (22 up_in_ONFH + DCN up_in_Control).
#   Null statistic used : fraction of random-set genes with logFC > 0 (up-fraction).
#   MISMATCH            : the null does not model the frozen direction pattern
#                         (23:1 up:down). Symmetric correction: assign the same
#                         direction pattern (22 up, 1 down; "down" position drawn
#                         uniformly at random) to each expression-matched random
#                         set and compute the identical match-rate statistic.
#   The ssGSEA disease-effect statistic (observed vs random sets) is already
#   definitionally symmetric (same set-level statistic on both) and is verified,
#   not recomputed.
#
# 3.2 Effect-size comparability
#   Raw ssGSEA disease effects are NOT comparable across signatures of different
#   size (Fibrocartilage 23 genes vs EMT 200 / ECM 111 / Fibroblast 88).
#   Added interpretable, size-robust references computed from the SAVED
#   per-sample scores: Cohen's d (disease vs control) for every signature, and
#   a null-standardized effect z for Fibrocartilage from the saved 1,000-set
#   null. Three significance types are reported separately.
#
# Output: tables/R2_extspec_stats_review.csv,
#         tables/R2_extspec_corrected_concordance.csv
# Seed 42. Existing R2_extspec_* files untouched.
# ==============================================================================
set.seed(42)
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages({library(GSVA); library(org.Hs.eg.db); library(AnnotationDbi); library(openxlsx)})
TAB <- "/mnt/results/phase1/tables"
GEO <- "/mnt/shared-workspace/shared/phase1/geo_data"
N_PERM <- 1000

cat("== [1] Rebuild expression matrices (deterministic, same code as original) ==\n")
expr74089_probe <- readRDS(file.path(GEO, "GSE74089_expr.rds"))
probe_annot <- readRDS(file.path(GEO, "GSE74089_probe_annotation.rds"))
annot <- probe_annot[probe_annot$probe %in% rownames(expr74089_probe) & !is.na(probe_annot$symbol), ]
expr_annot <- expr74089_probe[annot$probe, , drop = FALSE]
annot$variance <- apply(expr_annot, 1, var)
annot <- annot[order(annot$symbol, -annot$variance), ]
annot_dedup <- annot[!duplicated(annot$symbol), ]
X74089 <- as.matrix(expr74089_probe[annot_dedup$probe, , drop = FALSE])
rownames(X74089) <- annot_dedup$symbol

tpm <- read.xlsx(file.path(GEO, "GSE345158_gene_tpm_exp.annotation.xlsx"))
tpm_mat <- as.matrix(apply(tpm[, -1, drop = FALSE], 2, as.numeric))
rownames(tpm_mat) <- tpm$gene_id
colnames(tpm_mat) <- colnames(tpm)[-1]
sym <- mapIds(org.Hs.eg.db, keys = rownames(tpm_mat), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
keep <- !is.na(sym); tpm_mat <- tpm_mat[keep, ]; sym <- sym[keep]
ltpm <- log2(tpm_mat + 1)
ord <- order(sym, -apply(ltpm, 1, var)); sym_ord <- sym[ord]; ltpm_ord <- ltpm[ord, ]
keep2 <- !duplicated(sym_ord)
X345158 <- ltpm_ord[keep2, ]; rownames(X345158) <- sym_ord[keep2]

cat("== [2] Regenerate 1,000 expression-matched random sets (seed 42) ==\n")
frozen <- read.csv(file.path(TAB, "frozen_signatures.csv"))
fibro <- frozen[frozen$signature_name == "Fibrocartilage", ]
sig_genes <- fibro$gene

make_random_sets <- function(X, sig_genes, n_perm = N_PERM, seed = 42) {
  set.seed(seed)
  genes <- rownames(X); gm <- rowMeans(X)
  bins <- cut(rank(gm), breaks = 10, labels = FALSE); names(bins) <- genes
  sig_in <- intersect(sig_genes, genes); sig_bins <- bins[sig_in]
  bin_counts <- table(sig_bins)
  sets <- vector("list", n_perm)
  for (i in seq_len(n_perm)) {
    sets[[i]] <- unlist(lapply(names(bin_counts), function(b) {
      cand <- setdiff(names(bins)[bins == as.integer(b)], names(sig_bins))
      if (length(cand) >= bin_counts[b]) sample(cand, bin_counts[b]) else
        sample(cand, bin_counts[b], replace = TRUE)
    }))
  }
  sets
}

de74089 <- read.csv(file.path(TAB, "validation_GSE74089_DE.csv"))
de345158 <- read.csv(file.path(TAB, "validation_GSE345158_DE.csv"))
lfc7 <- de74089$logFC; names(lfc7) <- de74089$gene
lfc3 <- de345158$logFC; names(lfc3) <- de345158$gene
null_saved <- read.csv(file.path(TAB, "R2_extspec_random_null_raw.csv"))

review <- list(); conc_out <- list()
for (cfg in list(list(cohort = "GSE74089", X = X74089, lfc = lfc7),
                 list(cohort = "GSE345158", X = X345158, lfc = lfc3))) {
  co <- cfg$cohort; lfc <- cfg$lfc
  sets <- make_random_sets(cfg$X, sig_genes)
  sig_present <- intersect(sig_genes, rownames(cfg$X))
  expected_up <- fibro$direction == "up_in_ONFH"; names(expected_up) <- fibro$gene
  expected_up <- expected_up[sig_present]
  n_up <- sum(expected_up); n_down <- length(expected_up) - n_up

  # original null statistic (up-fraction) — verify against saved raw table
  conc_orig <- sapply(sets, function(g) { gg <- intersect(g, names(lfc)); mean(lfc[gg] > 0) })
  sv <- null_saved$concordance_rate[null_saved$cohort == co]
  regen_ok <- isTRUE(all.equal(round(conc_orig, 6), round(sv, 6)))
  cat(sprintf("[%s] random-set regeneration reproduces saved null: %s (max abs diff %.2e)\n",
              co, regen_ok, max(abs(conc_orig - sv))))

  # observed concordance (direction-match rate)
  conc_obs <- mean((lfc[sig_present] > 0) == expected_up)

  # CORRECTED symmetric null: same direction pattern (n_up up, n_down down),
  # "down" position drawn uniformly at random per random set
  set.seed(42 + 1)
  conc_sym <- sapply(sets, function(g) {
    gg <- intersect(g, names(lfc))
    if (length(gg) != length(sig_present)) return(NA)
    down_pos <- sample.int(length(gg), n_down)
    eu <- rep(TRUE, length(gg)); eu[down_pos] <- FALSE
    mean((lfc[gg] > 0) == eu)
  })
  p_orig <- (sum(conc_orig >= conc_obs, na.rm = TRUE) + 1) / (sum(!is.na(conc_orig)) + 1)
  p_sym  <- (sum(conc_sym >= conc_obs, na.rm = TRUE) + 1) / (sum(!is.na(conc_sym)) + 1)

  # second symmetric reference: direction-label permutation within the observed signature
  set.seed(42 + 2)
  perm <- replicate(N_PERM, {
    eu <- rep(TRUE, length(sig_present))
    eu[sample.int(length(sig_present), n_down)] <- FALSE
    mean((lfc[sig_present] > 0) == eu)
  })
  p_perm <- (sum(perm >= conc_obs) + 1) / (N_PERM + 1)

  # verify disease-effect empirical p (definitionally symmetric; recompute from saved null)
  tab1 <- read.csv(file.path(TAB, "R2_extspec_random_results.csv"))
  t1 <- tab1[tab1$cohort == co, ]
  nd <- null_saved$disease_effect[null_saved$cohort == co]
  p_eff_recomp <- (sum(abs(nd) >= abs(t1$observed_effect)) + 1) / (length(nd) + 1)
  eff_sym_ok <- isTRUE(all.equal(p_eff_recomp, t1$empirical_p_effect))
  z_eff <- (t1$observed_effect - t1$null_mean) / t1$null_sd

  review[[length(review) + 1]] <- data.frame(
    cohort = co, item = "directional_concordance",
    original = sprintf("observed=%.3f (%d/%d); null=up-fraction mean=%.3f; empirical p=%.4f",
                       conc_obs, round(conc_obs * length(sig_present)), length(sig_present),
                       mean(conc_orig, na.rm = TRUE), p_orig),
    corrected = sprintf("symmetric null (%d-up/%d-down pattern on matched random sets) mean=%.3f; p=%.4f; direction-permutation null p=%.4f",
                        n_up, n_down, mean(conc_sym, na.rm = TRUE), p_sym, p_perm),
    verdict = sprintf("regeneration faithful: %s; %s", regen_ok,
                      ifelse(p_sym <= 0.05, "concordance exceeds symmetric null (conclusion unchanged)",
                             "concordance NOT significant vs symmetric null (conclusion changes)")))
  review[[length(review) + 1]] <- data.frame(
    cohort = co, item = "ssGSEA_disease_effect_vs_random_null",
    original = sprintf("observed=%.3f; null mean=%.3f sd=%.3f; empirical p=%.4f (two-sided, add-one); percentile=%.1f",
                       t1$observed_effect, t1$null_mean, t1$null_sd, t1$empirical_p_effect, t1$percentile_effect),
    corrected = sprintf("definition already symmetric (same set-level statistic on observed and random sets); recomputed p=%.4f, identical to saved=%s; standardized z=%.2f",
                        p_eff_recomp, eff_sym_ok, z_eff),
    verdict = "no correction required")

  # 3.2 standardized effect for Fibrocartilage (saved null) + Cohen's d (saved scores)
  ps <- read.csv(file.path(TAB, "R2_extspec_ssgsea_scores.csv"))
  ps_c <- if (co == "GSE74089") ps[grepl("^GSM", ps$sample), ] else ps[!grepl("^GSM", ps$sample), ]
  is_d <- ps_c$group == "disease"; is_c <- ps_c$group == "control"
  sigs <- c("Fibrocartilage","Hallmark_EMT","Generic_ECM","Generic_Fibroblast")
  d_out <- sapply(sigs, function(s) {
    x1 <- ps_c[is_d, s]; x0 <- ps_c[is_c, s]
    (mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2)
  })
  eff_raw <- sapply(sigs, function(s) mean(ps_c[is_d, s]) - mean(ps_c[is_c, s]))
  review[[length(review) + 1]] <- data.frame(
    cohort = co, item = "effect_size_comparability",
    original = sprintf("raw ssGSEA disease effects: Fibrocartilage(23g)=%.3f, EMT(200g)=%.3f, ECM(111g)=%.3f, Fibroblast(88g)=%.3f",
                       eff_raw[["Fibrocartilage"]], eff_raw[["Hallmark_EMT"]],
                       eff_raw[["Generic_ECM"]], eff_raw[["Generic_Fibroblast"]]),
    corrected = sprintf("Cohen's d (size-robust): Fibrocartilage=%.2f, EMT=%.2f, ECM=%.2f, Fibroblast=%.2f; Fibrocartilage null-standardized z=%.2f",
                        d_out[["Fibrocartilage"]], d_out[["Hallmark_EMT"]], d_out[["Generic_ECM"]],
                        d_out[["Generic_Fibroblast"]], z_eff),
    verdict = "raw effects not comparable across set sizes; use d / z and the three significance types separately")
  conc_out[[length(conc_out) + 1]] <- data.frame(
    cohort = co, observed_match_rate = conc_obs, n_genes = length(sig_present),
    n_up_expected = n_up, n_down_expected = n_down,
    null_upfraction_mean = mean(conc_orig, na.rm = TRUE), p_original = p_orig,
    null_symmetric_mean = mean(conc_sym, na.rm = TRUE), p_corrected_symmetric = p_sym,
    p_corrected_direction_permutation = p_perm, regeneration_faithful = regen_ok)
}
out <- do.call(rbind, review)
print(out)
write.csv(out, file.path(TAB, "R2_extspec_stats_review.csv"), row.names = FALSE)
cout <- do.call(rbind, conc_out)
write.csv(cout, file.path(TAB, "R2_extspec_corrected_concordance.csv"), row.names = FALSE)
cat("\nDONE extspec stats review\n")
