#!/usr/bin/env R
# =============================================================================
# R2 Task 1 (revision round 2): Doublet sensitivity analysis with
# cartilage-lineage comparator
# =============================================================================
# Reviewer 2: compare doublet fraction between fibrocartilage-dominant cells
# and OTHER CARTILAGE-LINEAGE cells (not "all other cells"), then repeat the
# donor-level fibrocartilage score comparison after removing predicted doublets.
#
# Comparator definitions (from original pipeline metadata):
#   A  Fibrocartilage-dominant : dominant_program == "Fibrocartilage"
#   B  Other cartilage-lineage : lineage_call == "cartilage_lineage"
#                                & dominant_program != "Fibrocartilage"
#   B2 Sensitivity (stricter)  : B excluding dominant_program == "Undetermined"
#   C  Previous "Other_cells"  : all non-fibrocartilage cells (reported for
#                                continuity with R2_dblsens_* results only)
#
# lineage_call provenance (phase 0.5, cluster-level): clusters with supportive
# chondrocyte markers and NO fibroblast/stromal supportive markers were called
# "cartilage_lineage"; clusters with BOTH chondrocyte and fibroblast/stromal
# markers were called "cartilage_fibrocartilage_candidate"
# (see /mnt/shared-workspace/shared/phase0_5/GSE311084_lineage.rds).
#
# CORRECTED (R2 integrity audit, 2026-10-15): the earlier comment claimed all
# cells are "phase-0 Chondrocyte calls". That was imprecise. Phase-0 per-sample
# preliminary annotations were heterogeneous (Fibroblast 38,585 / Chondrocyte
# 18,765 / Macrophage 126 of 57,476 QC'd cells). All 57,310 RETAINED cells are
# Chondrocyte calls from the PHASE-0.5 integrated re-annotation (cluster-level,
# 9-lineage marker panel); the 166 cells in the phase-0.5 Macrophage cluster
# (ONFH-1: 7, ONFH-2: 159) were removed before the final object. See
# R2_cell_identity_audit.md and tables/R2_cell_annotation_provenance.csv.
# dominant_program "Undetermined" reflects absence of a dominant remodeling
# program, NOT lineage uncertainty.
#
# Seed 42. Donor = biological unit. Existing R2_dblsens_* results preserved.
# =============================================================================

set.seed(42)
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages({
  library(Seurat); library(ggplot2)
})
OUT_TAB <- "/mnt/results/phase1/tables"
OUT_FIG <- "/mnt/results/phase1/figures"
OUT_LOG <- "/mnt/results/phase1/logs"
dir.create(OUT_LOG, showWarnings = FALSE, recursive = TRUE)

cat("=== R2 Task 1 (round 2): doublet cartilage-lineage comparator ===\n")
cat("R version:", R.version.string, "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("Run started:", format(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load object with doublet calls (no re-preprocessing)
# -----------------------------------------------------------------------------
seu <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored_doublet.rds")
meta <- seu@meta.data
cat("Cells:", nrow(meta), "\n")
cat("scDblFinder classes:", paste(names(table(meta$scDblFinder_class)), collapse = ", "), "\n\n")

meta$is_doublet <- meta$scDblFinder_class == "doublet"
meta$is_fibro_dominant <- meta$dominant_program == "Fibrocartilage"

# Group definitions
meta$grp <- NA
meta$grp[meta$is_fibro_dominant] <- "A_Fibrocartilage_dominant"
meta$grp[!meta$is_fibro_dominant & meta$lineage_call == "cartilage_lineage"] <- "B_Other_cartilage_lineage"
meta$grp[!meta$is_fibro_dominant & meta$lineage_call == "cartilage_fibrocartilage_candidate"] <- "C_Other_fibrocartilage_candidate"

cat("== Group sizes ==\n")
print(table(meta$grp, useNA = "ifany"))
cat("\n== Group B program composition (cartilage-lineage, non-fibrocartilage) ==\n")
print(table(meta$dominant_program[meta$grp == "B_Other_cartilage_lineage"]))

# Sensitivity sub-comparator: B excluding program-Undetermined
meta$grp_sens <- meta$grp
meta$grp_sens[meta$grp == "B_Other_cartilage_lineage" &
                meta$dominant_program != "Undetermined"] <- "B2_B_without_Undetermined"

donors <- c("Ctrl-1", "Ctrl-2", "Ctrl-3", "ONFH-1", "ONFH-2", "ONFH-3")
cond_map <- setNames(c("Control", "Control", "Control", "ONFH", "ONFH", "ONFH"), donors)
meta$donor_id <- factor(meta$donor_id, levels = donors)

# -----------------------------------------------------------------------------
# 2. Doublet rates: per donor x group + totals
# -----------------------------------------------------------------------------
rate_fun <- function(m, grp_col, grp_name, d) {
  sel <- m[[grp_col]] == grp_name & m$donor_id == d
  n <- sum(sel); nd <- sum(sel & m$is_doublet)
  data.frame(group_definition = grp_name, donor_id = d, condition = cond_map[[d]],
             n_cells = n, n_doublets = nd,
             doublet_pct = if (n > 0) 100 * nd / n else NA_real_,
             note = if (n == 0) "no cells in this group for this donor" else "")
}

grp_defs <- list(
  c("grp", "A_Fibrocartilage_dominant"),
  c("grp", "B_Other_cartilage_lineage"),
  c("grp_sens", "B2_B_without_Undetermined"),
  c("grp", "C_Other_fibrocartilage_candidate")
)

per_donor <- do.call(rbind, lapply(grp_defs, function(gd) {
  do.call(rbind, lapply(donors, function(d) rate_fun(meta, gd[1], gd[2], d)))
}))
cat("\n== Per-donor doublet rates ==\n"); print(per_donor, digits = 3)

# Group totals
totals <- do.call(rbind, lapply(grp_defs, function(gd) {
  sel <- meta[[gd[1]]] == gd[2]
  n <- sum(sel); nd <- sum(sel & meta$is_doublet)
  data.frame(group_definition = gd[2], n_cells = n, n_doublets = nd,
             doublet_pct = if (n > 0) 100 * nd / n else NA_real_,
             note = if (n == 0) "no cells" else "")
}))
cat("\n== Group totals ==\n"); print(totals, digits = 3)

write.csv(per_donor, file.path(OUT_TAB, "R2_dblsens_cartilage_lineage_per_donor.csv"), row.names = FALSE)

# Combined rates table (totals for all definitions + previous Other_cells reference)
sel_prev <- !meta$is_fibro_dominant  # previous "Other_cells" definition
prev_tot <- data.frame(group_definition = "REF_previous_Other_cells_all_non_fibro",
                       n_cells = sum(sel_prev),
                       n_doublets = sum(sel_prev & meta$is_doublet),
                       doublet_pct = 100 * sum(sel_prev & meta$is_doublet) / sum(sel_prev),
                       note = "reference: previous R2_dblsens definition (all non-fibrocartilage cells)")
rates <- rbind(totals, prev_tot)
write.csv(rates, file.path(OUT_TAB, "R2_dblsens_cartilage_lineage_rates.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Statistics: Fisher exact (descriptive, cell-level) + paired donor Wilcoxon
# -----------------------------------------------------------------------------
tests <- data.frame()
# A vs B (primary, reviewer-requested comparison)
mk2x2 <- function(m, g1, g2, grp_col = "grp") {
  sel1 <- m[[grp_col]] == g1; sel2 <- m[[grp_col]] == g2
  matrix(c(sum(sel1 & m$is_doublet), sum(sel1 & !m$is_doublet),
           sum(sel2 & m$is_doublet), sum(sel2 & !m$is_doublet)),
         nrow = 2, dimnames = list(c("doublet", "singlet"), c(g1, g2)))
}
tab_AB <- mk2x2(meta, "A_Fibrocartilage_dominant", "B_Other_cartilage_lineage")
f_AB <- fisher.test(tab_AB)
cat("\n== Fisher exact A vs B (cell-level, descriptive) ==\n"); print(tab_AB)
cat(sprintf("p = %.3g, OR = %.4f\n", f_AB$p.value, f_AB$estimate))

# A vs B2 (sensitivity)
tab_AB2 <- mk2x2(meta, "A_Fibrocartilage_dominant", "B2_B_without_Undetermined", "grp_sens")
f_AB2 <- fisher.test(tab_AB2)
cat(sprintf("A vs B2 (excl. Undetermined): p = %.3g, OR = %.4f\n", f_AB2$p.value, f_AB2$estimate))

# Paired donor-level Wilcoxon (A vs B), with per-donor paired differences
piv <- reshape(per_donor[per_donor$group_definition %in%
                           c("A_Fibrocartilage_dominant", "B_Other_cartilage_lineage"), ],
               idvar = "donor_id", timevar = "group_definition",
               drop = c("condition", "note"),
               direction = "wide")
colnames(piv) <- gsub("doublet_pct\\.", "pct_", colnames(piv))
piv$paired_diff_pct <- piv$pct_A_Fibrocartilage_dominant - piv$pct_B_Other_cartilage_lineage
cat("\n== Per-donor paired differences (A pct - B pct) ==\n")
print(piv[, c("donor_id", "n_cells.A_Fibrocartilage_dominant", "pct_A_Fibrocartilage_dominant",
              "n_cells.B_Other_cartilage_lineage", "pct_B_Other_cartilage_lineage",
              "paired_diff_pct")], digits = 3)

wt <- wilcox.test(piv$pct_A_Fibrocartilage_dominant, piv$pct_B_Other_cartilage_lineage,
                  paired = TRUE, na.rm = TRUE)
cat(sprintf("\nPaired Wilcoxon (n=%d complete pairs): p = %.4f\n",
            sum(complete.cases(piv[, c("pct_A_Fibrocartilage_dominant",
                                       "pct_B_Other_cartilage_lineage")])), wt$p.value))

tests <- data.frame(
  comparison = c("A_vs_B_fisher_cell_level", "A_vs_B_fisher_OR",
                 "A_vs_B2exclUndetermined_fisher_cell_level", "A_vs_B2exclUndetermined_fisher_OR",
                 "A_vs_B_paired_wilcoxon_donor_level", "n_complete_pairs"),
  value = c(f_AB$p.value, f_AB$estimate, f_AB2$p.value, f_AB2$estimate, wt$p.value,
            sum(complete.cases(piv[, c("pct_A_Fibrocartilage_dominant",
                                       "pct_B_Other_cartilage_lineage")]))),
  note = c("descriptive cell-level test; donor is the biological unit",
           "odds of doublet in A vs B", "sensitivity comparator",
           "odds of doublet in A vs B2", "paired by donor, n=6",
           "donors with cells in both groups"))
write.csv(tests, file.path(OUT_TAB, "R2_dblsens_cartilage_lineage_fraction_tests.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Score sensitivity: donor-level fibrocartilage score before/after doublet removal
#    Primary: MS_2 (AddModuleScore; identical to manuscript donor_module_scores.csv)
#    Secondary: FibrocartilageUC_ (UCell)
# -----------------------------------------------------------------------------
score_col <- "MS_2"; score_col_uc <- "FibrocartilageUC_"
stopifnot(score_col %in% colnames(meta), score_col_uc %in% colnames(meta))

donor_stats <- data.frame(
  donor_id = donors, condition = cond_map[donors],
  ms2_all = sapply(donors, function(d) mean(meta[[score_col]][meta$donor_id == d])),
  ms2_singlet = sapply(donors, function(d) mean(meta[[score_col]][meta$donor_id == d & !meta$is_doublet])),
  ucell_all = sapply(donors, function(d) mean(meta[[score_col_uc]][meta$donor_id == d])),
  ucell_singlet = sapply(donors, function(d) mean(meta[[score_col_uc]][meta$donor_id == d & !meta$is_doublet])),
  n_all = sapply(donors, function(d) sum(meta$donor_id == d)),
  n_singlet = sapply(donors, function(d) sum(meta$donor_id == d & !meta$is_doublet))
)
cat("\n== Donor-level fibrocartilage scores (before/after doublet removal) ==\n")
print(donor_stats, digits = 3)

compare_groups <- function(score_vec, cond_vec) {
  x <- score_vec[cond_vec == "ONFH"]; y <- score_vec[cond_vec == "Control"]
  tt <- t.test(x, y, var.equal = FALSE)
  nx <- length(x); ny <- length(y)
  s <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  d <- (mean(x) - mean(y)) / s
  all_splits <- combn(seq_along(score_vec), 3)  # exact: C(6,3)=20 splits
  perm_diffs <- apply(all_splits, 2, function(idx) mean(score_vec[idx]) - mean(score_vec[-idx]))
  obs <- mean(x) - mean(y)
  list(effect = obs, welch_t = unname(tt$statistic), welch_p = tt$p.value,
       cohens_d = d, p_perm_one = mean(perm_diffs >= obs),
       p_perm_two = mean(abs(perm_diffs) >= abs(obs)))
}
set.seed(42)
mc_perm <- function(score_vec, cond_vec, n = 100000) {
  obs <- mean(score_vec[cond_vec == "ONFH"]) - mean(score_vec[cond_vec == "Control"])
  perms <- replicate(n, {
    idx <- sample(seq_along(score_vec), 3)
    mean(score_vec[idx]) - mean(score_vec[-idx])
  })
  (sum(abs(perms) >= abs(obs)) + 1) / (n + 1)
}

r1 <- compare_groups(donor_stats$ms2_all, donor_stats$condition)
r2 <- compare_groups(donor_stats$ms2_singlet, donor_stats$condition)
r3 <- compare_groups(donor_stats$ucell_all, donor_stats$condition)
r4 <- compare_groups(donor_stats$ucell_singlet, donor_stats$condition)

sens <- data.frame(
  score_type = rep(c("AddModuleScore_MS2_primary", "UCell_secondary"), each = 2),
  analysis = rep(c("before_doublet_removal", "after_doublet_removal"), 2),
  control_mean = c(mean(donor_stats$ms2_all[donor_stats$condition == "Control"]),
                   mean(donor_stats$ms2_singlet[donor_stats$condition == "Control"]),
                   mean(donor_stats$ucell_all[donor_stats$condition == "Control"]),
                   mean(donor_stats$ucell_singlet[donor_stats$condition == "Control"])),
  control_sd = c(sd(donor_stats$ms2_all[donor_stats$condition == "Control"]),
                 sd(donor_stats$ms2_singlet[donor_stats$condition == "Control"]),
                 sd(donor_stats$ucell_all[donor_stats$condition == "Control"]),
                 sd(donor_stats$ucell_singlet[donor_stats$condition == "Control"])),
  onfh_mean = c(mean(donor_stats$ms2_all[donor_stats$condition == "ONFH"]),
                mean(donor_stats$ms2_singlet[donor_stats$condition == "ONFH"]),
                mean(donor_stats$ucell_all[donor_stats$condition == "ONFH"]),
                mean(donor_stats$ucell_singlet[donor_stats$condition == "ONFH"])),
  onfh_sd = c(sd(donor_stats$ms2_all[donor_stats$condition == "ONFH"]),
              sd(donor_stats$ms2_singlet[donor_stats$condition == "ONFH"]),
              sd(donor_stats$ucell_all[donor_stats$condition == "ONFH"]),
              sd(donor_stats$ucell_singlet[donor_stats$condition == "ONFH"])),
  effect_direction = ifelse(c(r1$effect, r2$effect, r3$effect, r4$effect) > 0,
                            "higher_in_ONFH", "higher_in_Control"),
  welch_t = c(r1$welch_t, r2$welch_t, r3$welch_t, r4$welch_t),
  welch_p = c(r1$welch_p, r2$welch_p, r3$welch_p, r4$welch_p),
  cohens_d = c(r1$cohens_d, r2$cohens_d, r3$cohens_d, r4$cohens_d),
  perm_p_exact_one_sided = c(r1$p_perm_one, r2$p_perm_one, r3$p_perm_one, r4$p_perm_one),
  perm_p_exact_two_sided = c(r1$p_perm_two, r2$p_perm_two, r3$p_perm_two, r4$p_perm_two),
  perm_p_montecarlo_two_sided = c(mc_perm(donor_stats$ms2_all, donor_stats$condition),
                                  mc_perm(donor_stats$ms2_singlet, donor_stats$condition),
                                  NA, NA),
  n_donors = 3, n_perm_exact = 20,
  note = "exact permutation enumerates all C(6,3)=20 donor-label splits; minimum two-sided p = 0.10 at n=3 vs 3"
)
cat("\n== Score sensitivity summary ==\n"); print(sens, digits = 4)
write.csv(sens, file.path(OUT_TAB, "R2_dblsens_cartilage_lineage_score_sensitivity.csv"), row.names = FALSE)
write.csv(donor_stats, file.path(OUT_TAB, "R2_dblsens_cartilage_lineage_donor_scores.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Figure: A = per-donor doublet fraction (A vs B); B = before/after scores
# -----------------------------------------------------------------------------
svg_ok <- requireNamespace("svglite", quietly = TRUE)
save_fig <- function(p, name, w = 9.5, h = 4) {
  if (svg_ok) svglite::svglite(file.path(OUT_FIG, paste0(name, ".svg")), width = w, height = h)
  else svg(file.path(OUT_FIG, paste0(name, ".svg")), width = w, height = h)
  print(p); dev.off()
  ggsave(file.path(OUT_FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

pd <- per_donor[per_donor$group_definition %in%
                  c("A_Fibrocartilage_dominant", "B_Other_cartilage_lineage"), ]
pd$donor_id <- factor(pd$donor_id, levels = donors)
pd$group_definition <- factor(pd$group_definition,
                              levels = c("A_Fibrocartilage_dominant", "B_Other_cartilage_lineage"),
                              labels = c("Fibrocartilage-dominant", "Other cartilage-lineage"))
pA <- ggplot(pd, aes(x = donor_id, y = doublet_pct, fill = group_definition)) +
  geom_col(position = position_dodge(0.75), width = 0.7, color = "grey25", linewidth = 0.2) +
  scale_fill_manual(values = c("Fibrocartilage-dominant" = "#0279EE",
                               "Other cartilage-lineage" = "#ECE9E2"), name = NULL) +
  labs(x = NULL, y = "Predicted doublets (%)", title = "A") +
  theme_classic(base_family = "Liberation Sans") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

ps <- donor_stats
ps$condition <- factor(ps$condition, levels = c("Control", "ONFH"))
pl <- rbind(
  data.frame(donor = ps$donor_id, condition = ps$condition,
             stage = "All cells", score = ps$ms2_all),
  data.frame(donor = ps$donor_id, condition = ps$condition,
             stage = "Doublets removed", score = ps$ms2_singlet))
pl$stage <- factor(pl$stage, levels = c("All cells", "Doublets removed"))
pB <- ggplot(pl, aes(x = stage, y = score, fill = condition, group = donor)) +
  geom_point(position = position_dodge(0.25), size = 2.4, color = "grey25",
             shape = 21, stroke = 0.4) +
  geom_line(position = position_dodge(0.25), color = "grey55", linewidth = 0.3) +
  stat_summary(aes(group = stage), fun = mean, geom = "crossbar",
               width = 0.4, linewidth = 0.4, color = "black") +
  scale_fill_manual(values = c(Control = "#9DBCD4", ONFH = "#FF9400"), name = NULL) +
  labs(x = NULL, y = "Donor-level Fibrocartilage\nmodule score (AddModuleScore)", title = "B") +
  theme_classic(base_family = "Liberation Sans") +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

library(patchwork)
save_fig(pA | pB, "R2_Fig_dblsens_cartilage_lineage", w = 9.5, h = 4.2)

writeLines(capture.output(sessionInfo()),
           file.path(OUT_TAB, "R2_dblsens_cartilage_lineage_sessionInfo.txt"))
cat("\nRun finished:", format(Sys.time()), "\n")
cat("=== R2 Task 1 (round 2) Complete ===\n")
