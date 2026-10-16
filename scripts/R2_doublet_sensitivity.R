#!/usr/bin/env R
# =============================================================================
# R2 Analysis 2: Doublet sensitivity analysis (GSE311084)
# =============================================================================
# Purpose: Assess whether the fibrocartilage phenotype is driven by predicted
# heterotypic doublets (scDblFinder calls from revision Task 3; NOT re-run).
# Part 1: doublet fraction, Fibrocartilage-dominant vs other cells
# Part 2: donor-level Fibrocartilage UCell score, before vs after doublet removal
# Seeds: 42 (permutation). No re-preprocessing; existing scored object reused.
# =============================================================================

set.seed(42)
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages({
  library(Seurat); library(ggplot2)
})
OUT_TAB <- "/mnt/results/phase1/tables"
OUT_FIG <- "/mnt/results/phase1/figures"
N_PERM <- 100000  # random label permutations (exhaustive 20 splits reported too)

cat("=== R2 Analysis 2: Doublet sensitivity ===\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load scored object with doublet calls (no re-preprocessing)
# -----------------------------------------------------------------------------
seu <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored_doublet.rds")
meta <- seu@meta.data
cat("Cells:", nrow(meta), "\n")
cat("scDblFinder classes:", paste(names(table(meta$scDblFinder_class)), collapse = ", "), "\n")
cat("dominant_program levels:", paste(levels(meta$dominant_program), collapse = ", "), "\n")
cat("lineage_call values:", paste(names(table(meta$lineage_call)), collapse = ", "), "\n\n")

score_col <- "MS_2"  # AddModuleScore Fibrocartilage: matches manuscript donor_module_scores.csv
score_col_ucell <- "FibrocartilageUC_"  # secondary (UCell) sensitivity check
stopifnot(score_col %in% colnames(meta), score_col_ucell %in% colnames(meta))

meta$is_doublet <- meta$scDblFinder_class == "doublet"
meta$is_fibro_dominant <- meta$dominant_program == "Fibrocartilage"

donors <- c("Ctrl-1", "Ctrl-2", "Ctrl-3", "ONFH-1", "ONFH-2", "ONFH-3")
cond_map <- setNames(c("Control", "Control", "Control", "ONFH", "ONFH", "ONFH"), donors)

# -----------------------------------------------------------------------------
# 2. Part 1: doublet fraction, Fibrocartilage-dominant vs other cells
# -----------------------------------------------------------------------------
grp <- ifelse(meta$is_fibro_dominant, "Fibrocartilage_dominant", "Other_cells")
tab <- table(grp, meta$is_doublet)
cat("--- Doublet fraction (cell-level counts) ---\n"); print(tab)

fisher_res <- fisher.test(tab)
cat(sprintf("Fisher exact p = %.3g (OR = %.3f)\n",
    fisher_res$p.value, fisher_res$estimate))

# per-donor paired comparison (donor = biological unit, n=6)
donor_frac <- data.frame(
  donor_id = unique(meta$donor_id),
  stringsAsFactors = FALSE
)
donor_frac$fibro_n <- sapply(donor_frac$donor_id,
  function(d) sum(meta$donor_id == d & meta$is_fibro_dominant))
donor_frac$fibro_dbl <- sapply(donor_frac$donor_id,
  function(d) sum(meta$donor_id == d & meta$is_fibro_dominant & meta$is_doublet))
donor_frac$other_n <- sapply(donor_frac$donor_id,
  function(d) sum(meta$donor_id == d & !meta$is_fibro_dominant))
donor_frac$other_dbl <- sapply(donor_frac$donor_id,
  function(d) sum(meta$donor_id == d & !meta$is_fibro_dominant & meta$is_doublet))
donor_frac$fibro_pct <- 100 * donor_frac$fibro_dbl / donor_frac$fibro_n
donor_frac$other_pct <- 100 * donor_frac$other_dbl / donor_frac$other_n
donor_frac <- donor_frac[match(donors, donor_frac$donor_id), ]  # align to fixed donor order
cat("\n--- Per-donor doublet fraction (%) ---\n"); print(donor_frac, digits = 3)

wt_res <- wilcox.test(donor_frac$fibro_pct, donor_frac$other_pct, paired = TRUE)
cat(sprintf("Paired Wilcoxon (n=6 donors) p = %.3f\n", wt_res$p.value))

# cell composition cross-tab for transparency (lineage_call x group)
cat("\n--- lineage_call x group cross-tab ---\n")
print(table(meta$lineage_call, grp))

# -----------------------------------------------------------------------------
# 3. Part 2: donor-level Fibrocartilage score before vs after doublet removal
#    Primary: MS_2 (AddModuleScore, same as manuscript donor_module_scores.csv)
#    Secondary: FibrocartilageUC_ (UCell)
# -----------------------------------------------------------------------------
donor_stats <- data.frame(
  donor_id = donors,
  condition = cond_map[donors],
  score_all = sapply(donors, function(d)
    mean(meta[[score_col]][meta$donor_id == d])),
  score_singlet = sapply(donors, function(d)
    mean(meta[[score_col]][meta$donor_id == d & !meta$is_doublet])),
  score_all_ucell = sapply(donors, function(d)
    mean(meta[[score_col_ucell]][meta$donor_id == d])),
  score_singlet_ucell = sapply(donors, function(d)
    mean(meta[[score_col_ucell]][meta$donor_id == d & !meta$is_doublet])),
  n_all = sapply(donors, function(d) sum(meta$donor_id == d)),
  n_singlet = sapply(donors, function(d) sum(meta$donor_id == d & !meta$is_doublet))
)
cat("\n--- Donor-level Fibrocartilage scores ---\n"); print(donor_stats, digits = 3)

compare_groups <- function(score_vec, cond_vec) {
  x <- score_vec[cond_vec == "ONFH"]; y <- score_vec[cond_vec == "Control"]
  tt <- t.test(x, y, var.equal = FALSE)
  nx <- length(x); ny <- length(y)
  s <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  d <- (mean(x) - mean(y)) / s
  # exact permutation: all C(6,3)=20 label splits
  all_splits <- combn(seq_along(score_vec), 3)
  perm_diffs <- apply(all_splits, 2, function(idx) {
    mean(score_vec[idx]) - mean(score_vec[-idx])
  })
  obs <- mean(x) - mean(y)
  p_one <- mean(perm_diffs >= obs)          # includes observed split
  p_two <- mean(abs(perm_diffs) >= abs(obs))
  list(t = unname(tt$statistic), p_welch = tt$p.value, cohens_d = d,
       effect = obs, p_perm_one = p_one, p_perm_two = p_two,
       n_perm_exact = ncol(all_splits))
}

res_before <- compare_groups(donor_stats$score_all, donor_stats$condition)
res_after <- compare_groups(donor_stats$score_singlet, donor_stats$condition)
# UCell secondary
res_before_uc <- compare_groups(donor_stats$score_all_ucell, donor_stats$condition)
res_after_uc <- compare_groups(donor_stats$score_singlet_ucell, donor_stats$condition)

# Monte-Carlo permutation (N_PERM random label permutations) as robustness check
set.seed(42)
mc_perm <- function(score_vec, cond_vec, n = N_PERM) {
  obs <- mean(score_vec[cond_vec == "ONFH"]) - mean(score_vec[cond_vec == "Control"])
  perms <- replicate(n, {
    idx <- sample(seq_along(score_vec), 3)
    mean(score_vec[idx]) - mean(score_vec[-idx])
  })
  (sum(abs(perms) >= abs(obs)) + 1) / (n + 1)
}
res_before$p_perm_mc <- mc_perm(donor_stats$score_all, donor_stats$condition)
res_after$p_perm_mc <- mc_perm(donor_stats$score_singlet, donor_stats$condition)

sens <- data.frame(
  score_type = rep(c("AddModuleScore_MS2_primary", "UCell_secondary"), each = 2),
  analysis = rep(c("Before_doublet_removal", "After_doublet_removal"), 2),
  control_mean = c(mean(donor_stats$score_all[donor_stats$condition == "Control"]),
                   mean(donor_stats$score_singlet[donor_stats$condition == "Control"]),
                   mean(donor_stats$score_all_ucell[donor_stats$condition == "Control"]),
                   mean(donor_stats$score_singlet_ucell[donor_stats$condition == "Control"])),
  control_sd = c(sd(donor_stats$score_all[donor_stats$condition == "Control"]),
                 sd(donor_stats$score_singlet[donor_stats$condition == "Control"]),
                 sd(donor_stats$score_all_ucell[donor_stats$condition == "Control"]),
                 sd(donor_stats$score_singlet_ucell[donor_stats$condition == "Control"])),
  onfh_mean = c(mean(donor_stats$score_all[donor_stats$condition == "ONFH"]),
                mean(donor_stats$score_singlet[donor_stats$condition == "ONFH"]),
                mean(donor_stats$score_all_ucell[donor_stats$condition == "ONFH"]),
                mean(donor_stats$score_singlet_ucell[donor_stats$condition == "ONFH"])),
  onfh_sd = c(sd(donor_stats$score_all[donor_stats$condition == "ONFH"]),
              sd(donor_stats$score_singlet[donor_stats$condition == "ONFH"]),
              sd(donor_stats$score_all_ucell[donor_stats$condition == "ONFH"]),
              sd(donor_stats$score_singlet_ucell[donor_stats$condition == "ONFH"])),
  effect_direction = ifelse(c(res_before$effect, res_after$effect,
                              res_before_uc$effect, res_after_uc$effect) > 0,
                            "higher_in_ONFH", "higher_in_Control"),
  welch_t = c(res_before$t, res_after$t, res_before_uc$t, res_after_uc$t),
  welch_p = c(res_before$p_welch, res_after$p_welch,
              res_before_uc$p_welch, res_after_uc$p_welch),
  cohens_d = c(res_before$cohens_d, res_after$cohens_d,
               res_before_uc$cohens_d, res_after_uc$cohens_d),
  perm_p_one_sided = c(res_before$p_perm_one, res_after$p_perm_one,
                       res_before_uc$p_perm_one, res_after_uc$p_perm_one),
  perm_p_two_sided = c(res_before$p_perm_two, res_after$p_perm_two,
                       res_before_uc$p_perm_two, res_after_uc$p_perm_two),
  perm_p_montecarlo = c(res_before$p_perm_mc, res_after$p_perm_mc, NA, NA)
)
cat("\n--- Sensitivity summary ---\n"); print(sens, digits = 4)

# -----------------------------------------------------------------------------
# 4. Tables
# -----------------------------------------------------------------------------
# Doublet statistics table
dbl_tab <- data.frame(
  group = c(rep("Fibrocartilage_dominant", 6), rep("Other_cells", 6)),
  donor_id = rep(donors, 2),
  condition = rep(cond_map[donors], 2),
  n_cells = c(donor_frac$fibro_n, donor_frac$other_n),
  n_doublets = c(donor_frac$fibro_dbl, donor_frac$other_dbl),
  doublet_pct = c(donor_frac$fibro_pct, donor_frac$other_pct)
)
tot <- data.frame(group = c("Fibrocartilage_dominant", "Other_cells"),
  donor_id = "ALL", condition = "ALL",
  n_cells = c(sum(meta$is_fibro_dominant), sum(!meta$is_fibro_dominant)),
  n_doublets = c(sum(meta$is_fibro_dominant & meta$is_doublet),
                 sum(!meta$is_fibro_dominant & meta$is_doublet)))
tot$doublet_pct <- 100 * tot$n_doublets / tot$n_cells
dbl_tab <- rbind(dbl_tab, tot)
write.csv(dbl_tab, file.path(OUT_TAB, "R2_dblsens_doublet_stats.csv"), row.names = FALSE)

stats_row <- data.frame(
  test = c("Fisher_exact_cell_level_p", "Fisher_OR",
           "Paired_Wilcoxon_donor_level_p", "n_donors"),
  value = c(fisher_res$p.value, fisher_res$estimate, wt_res$p.value, 6)
)
write.csv(stats_row, file.path(OUT_TAB, "R2_dblsens_fraction_tests.csv"), row.names = FALSE)

write.csv(donor_stats, file.path(OUT_TAB, "R2_dblsens_donor_scores.csv"), row.names = FALSE)
write.csv(sens, file.path(OUT_TAB, "R2_dblsens_sensitivity_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Figure: doublet sensitivity (2 panels)
# -----------------------------------------------------------------------------
svg_ok <- requireNamespace("svglite", quietly = TRUE)
save_fig <- function(p, name, w = 9, h = 3.8) {
  if (svg_ok) svglite::svglite(file.path(OUT_FIG, paste0(name, ".svg")), width = w, height = h)
  else svg(file.path(OUT_FIG, paste0(name, ".svg")), width = w, height = h)
  print(p); dev.off()
  ggsave(file.path(OUT_FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

plot_dbl <- dbl_tab[dbl_tab$donor_id != "ALL", ]
plot_dbl$donor_id <- factor(plot_dbl$donor_id, levels = donors)
plot_dbl$group <- factor(plot_dbl$group, levels = c("Fibrocartilage_dominant", "Other_cells"))
pA <- ggplot(plot_dbl, aes(x = donor_id, y = doublet_pct, fill = group)) +
  geom_col(position = position_dodge(0.75), width = 0.7, color = "grey25", linewidth = 0.2) +
  scale_fill_manual(values = c(Fibrocartilage_dominant = "#0279EE", Other_cells = "#ECE9E2"),
                    name = NULL) +
  labs(x = NULL, y = "Predicted doublets (%)", title = "A") +
  theme_classic(base_family = "Liberation Sans") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

plot_sens <- donor_stats
plot_sens$condition <- factor(plot_sens$condition, levels = c("Control", "ONFH"))
plot_long <- rbind(
  data.frame(donor = plot_sens$donor_id, condition = plot_sens$condition,
             stage = "All cells", score = plot_sens$score_all),
  data.frame(donor = plot_sens$donor_id, condition = plot_sens$condition,
             stage = "Doublets removed", score = plot_sens$score_singlet))
plot_long$stage <- factor(plot_long$stage, levels = c("All cells", "Doublets removed"))
pB <- ggplot(plot_long, aes(x = stage, y = score, fill = condition, group = donor)) +
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
fig <- pA | pB
save_fig(fig, "R2_Fig_dblsens", w = 9.5, h = 4)

writeLines(capture.output(sessionInfo()),
           file.path(OUT_TAB, "R2_dblsens_sessionInfo.txt"))
cat("\n=== R2 Analysis 2 Complete ===\n")
