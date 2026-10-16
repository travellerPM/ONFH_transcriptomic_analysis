# ==============================================================================
# R2_task4_meta_analysis_audit.R
# Reviewer-2 integrity audit, Task 4: audit of the random-effects synthesis
# (R1_meta_analysis.R) and a more prudent descriptive synthesis.
#
# Audit findings (documented in R2_final_integrity_audit.md):
#   a) Effect-size metric = raw mean difference of ssGSEA scores (disease -
#      control). ssGSEA scores are within-sample rank-normalized (comparable
#      [0,1] scale across platforms), so a raw MD is defensible, but score SDs
#      differ between cohorts -> standardized effect (Hedges g) added as a
#      comparability check.
#   b) SE = |diff / t_stat| — algebraically exact recovery of the SE of the
#      mean difference from the per-cohort t-test. No issue.
#   c) DerSimonian-Laird estimator via metafor::rma(method="DL") — correctly
#      called; tau2 and I2 come from the same fit.
#   d) With k = 2 studies (df = 1), tau2/I2/Q are extremely unstable; I2 values
#      near 95-96% here are driven by single-df Q and must not be interpreted
#      substantively.
#   e) Prudent descriptive synthesis added: per-cohort effects + 95% CI,
#      direction concordance, fixed-effect vs random-effects comparison, and
#      Hedges-g SMD sensitivity. Inferential claims limited accordingly.
#
# Output: tables/R1_meta_audit_results.csv
# Seed 42. R1_meta_analysis_results.csv untouched (original preserved).
# ==============================================================================
set.seed(42)
suppressMessages(library(metafor))
TAB <- "/mnt/results/phase1/tables"

v7 <- read.csv(file.path(TAB, "validation_GSE74089_signature_scores.csv"))
v3 <- read.csv(file.path(TAB, "validation_GSE345158_signature_scores.csv"))
v7$dataset <- "GSE74089"; v3$dataset <- "GSE345158"
names(v3)[names(v3) == "sonfh_mean"] <- "nfh_mean"
v7$n_d <- 4; v7$n_c <- 4; v3$n_d <- 5; v3$n_c <- 4

comb <- rbind(
  v7[, c("signature","dataset","diff","t_stat","p_value","n_d","n_c")],
  v3[, c("signature","dataset","diff","t_stat","p_value","n_d","n_c")])
names(comb)[names(comb) == "diff"] <- "effect_size"
comb$se <- abs(comb$effect_size / comb$t_stat)

# Hedges g per cohort (standardized mean difference, small-sample corrected)
hedges_g <- function(diff, sd1, sd0, n1, n0) {
  s <- sqrt(((n1 - 1) * sd1^2 + (n0 - 1) * sd0^2) / (n1 + n0 - 2))
  d <- diff / s
  J <- 1 - 3 / (4 * (n1 + n0) - 9)
  d * J
}
# per-cohort SDs of ssGSEA scores recovered from saved per-sample scores
ps <- read.csv(file.path(TAB, "R2_extspec_ssgsea_scores.csv"))
ps7 <- ps[grepl("^GSM", ps$sample), ]; ps3 <- ps[!grepl("^GSM", ps$sample), ]
sig_names <- unique(comb$signature)

out <- data.frame()
for (sig in sig_names) {
  sub <- comb[comb$signature == sig, ]
  if (nrow(sub) < 2) next
  # reproduce original DL fit
  fit_dl <- rma(yi = effect_size, sei = se, data = sub, method = "DL", test = "z")
  fit_fe <- rma(yi = effect_size, sei = se, data = sub, method = "FE")
  # SMD (Hedges g) per cohort from per-sample ssGSEA scores
  # (available only for signatures in the saved per-sample score table:
  #  Fibrocartilage + 3 comparators; other signatures -> NA, MD synthesis only)
  g_vec <- c(NA_real_, NA_real_); g_se <- c(NA_real_, NA_real_)
  if (sig %in% colnames(ps7)) {
    for (ds in c("GSE74089", "GSE345158")) {
      psc <- if (ds == "GSE74089") ps7 else ps3
      x1 <- psc[psc$group == "disease", sig]; x0 <- psc[psc$group == "control", sig]
      n1 <- length(x1); n0 <- length(x0)
      g <- hedges_g(mean(x1) - mean(x0), sd(x1), sd(x0), n1, n0)
      se_g <- sqrt((n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0 - 2)))
      g_vec[if (ds == "GSE74089") 1 else 2] <- g
      g_se[if (ds == "GSE74089") 1 else 2] <- se_g
    }
  }
  smd_avail <- !is.na(g_vec[1]) && !is.na(g_vec[2])
  fit_smd_dl <- if (smd_avail) rma(yi = g_vec, sei = g_se, method = "DL", test = "z") else NULL
  dir_agree <- sign(sub$effect_size[1]) == sign(sub$effect_size[2])
  out <- rbind(out, data.frame(
    signature = sig,
    md_74089 = round(sub$effect_size[sub$dataset == "GSE74089"], 4),
    md_345158 = round(sub$effect_size[sub$dataset == "GSE345158"], 4),
    direction_concordant = dir_agree,
    dl_combined_md = round(fit_dl$beta[1], 4),
    dl_ci = paste0("[", round(fit_dl$ci.lb, 4), ", ", round(fit_dl$ci.ub, 4), "]"),
    dl_p = signif(fit_dl$pval, 3),
    fe_combined_md = round(fit_fe$beta[1], 4),
    fe_p = signif(fit_fe$pval, 3),
    Q_p = signif(fit_dl$QEp, 3),
    I2_pct = round(fit_dl$I2, 1),
    tau2 = signif(fit_dl$tau2, 3),
    hedges_g_74089 = if (smd_avail) round(g_vec[1], 3) else NA,
    hedges_g_345158 = if (smd_avail) round(g_vec[2], 3) else NA,
    smd_dl_combined_g = if (smd_avail) round(fit_smd_dl$beta[1], 3) else NA,
    smd_dl_ci = if (smd_avail) paste0("[", round(fit_smd_dl$ci.lb, 3), ", ", round(fit_smd_dl$ci.ub, 3), "]") else NA,
    smd_dl_p = if (smd_avail) signif(fit_smd_dl$pval, 3) else NA))
}
print(out, digits = 3)
write.csv(out, file.path(TAB, "R1_meta_audit_results.csv"), row.names = FALSE)

# verify reproduction of the original DL numbers
orig <- read.csv(file.path(TAB, "R1_meta_analysis_results.csv"))
m <- merge(orig[, c("signature","combined_effect")], out[, c("signature","dl_combined_md")],
           by = "signature")
cat("\nMax |original - recomputed| DL combined effect:",
    max(abs(m$combined_effect - m$dl_combined_md)), "\n")
cat("DONE meta audit\n")
