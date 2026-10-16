# ==============================================================================
# R2_task5_figure2_scoring_audit.R
# Reviewer-2 integrity audit, Task 5: Figure 2 scoring-system audit and
# re-output of consistent score figures.
#
# Audit findings:
#   - Fig2_module_scores.svg / Fig2b_umap_feature_scores.svg: CELL-level UCell
#     boxplots / UMAP feature plots (title and non-negative 0-1 scale confirm
#     UCell). Internally correctly labeled.
#   - Results Section 3.2 quotes 0.556 vs 1.830 = DONOR-level AddModuleScore
#     (MS_2) group means (donor_module_scores.csv). -> Algorithm (UCell vs
#     AddModuleScore) AND level (cell vs donor) differ between figure and text.
#   - R2_Fig_donor_scores*.svg: plots donor_module_scores.csv (AddModuleScore)
#     but labels the y axis "Donor mean UCell score" -> WRONG LABEL (data match
#     the text). Corrected versions emitted below; originals preserved.
#
# Re-outputs (originals untouched):
#   figures/R2_Fig_donor_scores_corrected.svg/.png      (label fix)
#   figures/R2_Fig_donor_scores_bar_corrected.svg/.png  (label fix)
#   figures/Fig2_module_scores_donor_AddModuleScore.svg/.png (matches text)
#   figures/Fig2_module_scores_donor_UCell.svg/.png          (UCell counterpart)
#   tables/R2_figure2_scoring_audit.csv
# Seed 42.
# ==============================================================================
set.seed(42)
suppressMessages({library(Seurat); library(ggplot2); library(patchwork)})

TAB <- "/mnt/results/phase1/tables"
FIG <- "/mnt/results/phase1/figures"

dm <- read.csv(file.path(TAB, "donor_module_scores.csv"))  # AddModuleScore donor means
dm$condition <- factor(dm$condition, levels = c("Control", "ONFH"))
dm$donor <- factor(dm$donor, levels = c("Ctrl-1","Ctrl-2","Ctrl-3","ONFH-1","ONFH-2","ONFH-3"))
prog_order <- c("Homeostatic","Fibrocartilage","Inflammatory","Pre_hypertrophic","Hypertrophic_osteogenic")
dm$program <- factor(dm$program, levels = prog_order)

theme_pub <- theme_bw(base_family = "Liberation Sans") +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "none")

save_fig <- function(p, name, w, h) {
  svglite::svglite(file.path(FIG, paste0(name, ".svg")), width = w, height = h)
  print(p); dev.off()
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}

# ---- 1. Corrected donor figure (label fix; data unchanged: AddModuleScore) ----
fc <- dm[dm$program == "Fibrocartilage", ]
tt <- t.test(mean_score ~ condition, data = fc)
fs <- aggregate(mean_score ~ condition, data = fc, FUN = mean)
d_cohen <- (fs$mean_score[fs$condition == "ONFH"] - fs$mean_score[fs$condition == "Control"]) /
  sqrt(sd(fc$mean_score[fc$condition == "Control"])^2 / 2 + sd(fc$mean_score[fc$condition == "ONFH"])^2 / 2)

p_a <- ggplot(fc, aes(x = condition, y = mean_score, fill = condition)) +
  geom_boxplot(width = 0.4, alpha = 0.3, outlier.shape = NA) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.8,
             position = position_jitter(width = 0.1, seed = 42)) +
  scale_fill_manual(values = c("Control" = "#0279EE", "ONFH" = "#FF9400")) +
  labs(title = "Fibrocartilage remodeling program",
       subtitle = sprintf("Donor-level mean AddModuleScore (MS_2), n=3 per group; Welch t-test p=%.4f, Cohen's d=%.2f",
                          tt$p.value, d_cohen),
       x = NULL, y = "Donor mean AddModuleScore") +
  theme_pub
p_b <- ggplot(dm, aes(x = condition, y = mean_score, fill = condition)) +
  geom_boxplot(width = 0.4, alpha = 0.3, outlier.shape = NA) +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.7,
             position = position_jitter(width = 0.1, seed = 42)) +
  facet_wrap(~program, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Control" = "#0279EE", "ONFH" = "#FF9400")) +
  labs(title = "All remodeling programs", x = NULL, y = "Donor mean AddModuleScore") +
  theme_pub
save_fig(p_a / p_b + plot_layout(heights = c(1, 1)), "R2_Fig_donor_scores_corrected", 14, 8)

p_bar <- ggplot(dm, aes(x = program, y = mean_score, fill = condition)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, alpha = 0.6) +
  geom_point(aes(group = condition), position = position_dodge(width = 0.6),
             size = 3, shape = 21, color = "black", stroke = 0.7) +
  scale_fill_manual(values = c("Control" = "#0279EE", "ONFH" = "#FF9400")) +
  labs(x = NULL, y = "Donor mean AddModuleScore",
       title = "Donor-level program scores (n=3 per group)") +
  theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
                    legend.position = "bottom", legend.title = element_blank())
save_fig(p_bar, "R2_Fig_donor_scores_bar_corrected", 10, 6)

# ---- 2. Donor-level AddModuleScore figure matching Results 3.2 (new Fig2 panel) ----
p_fc2 <- ggplot(fc, aes(x = condition, y = mean_score, fill = condition)) +
  geom_boxplot(width = 0.4, alpha = 0.3, outlier.shape = NA) +
  geom_point(size = 4.5, shape = 21, color = "black", stroke = 0.8,
             position = position_jitter(width = 0.08, seed = 42)) +
  geom_segment(data = data.frame(x1 = 1, x2 = 2, y = 2.35),
               aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE,
               linewidth = 0.4) +
  annotate("text", x = 1.5, y = 2.42, label = sprintf("p = %.4f", tt$p.value), size = 3.2) +
  scale_fill_manual(values = c("Control" = "#0279EE", "ONFH" = "#FF9400")) +
  scale_y_continuous(limits = c(0, 2.55)) +
  labs(title = "Fibrocartilage program, donor-level AddModuleScore",
       subtitle = "Donors as biological units (n=3 vs 3); error bars omitted, points are individual donors",
       x = NULL, y = "Donor mean AddModuleScore (MS_2)") +
  theme_pub
save_fig(p_fc2, "Fig2_module_scores_donor_AddModuleScore", 6, 5)

# ---- 3. Donor-level UCell counterpart (from scored object) ----
sc <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored.rds")
uc_cols <- c("HomeostaticUC_","FibrocartilageUC_","InflammatoryUC_",
             "Pre_hypertrophicUC_","Hypertrophic_osteogenicUC_")
names(uc_cols) <- prog_order
du <- do.call(rbind, lapply(prog_order, function(pr) {
  v <- sc@meta.data[[uc_cols[[pr]]]]
  m <- tapply(v, sc$donor_id, mean)
  cond <- tapply(as.character(sc$condition), sc$donor_id, function(x) unique(x)[1])
  data.frame(donor = names(m), condition = as.character(cond), program = pr,
             mean_score = as.numeric(m))
}))
du$condition <- factor(du$condition, levels = c("Control","ONFH"))
du$program <- factor(du$program, levels = prog_order)
tt_u <- t.test(mean_score ~ condition, data = du[du$program == "Fibrocartilage", ])
p_uc <- ggplot(du[du$program == "Fibrocartilage", ],
               aes(x = condition, y = mean_score, fill = condition)) +
  geom_boxplot(width = 0.4, alpha = 0.3, outlier.shape = NA) +
  geom_point(size = 4.5, shape = 21, color = "black", stroke = 0.8,
             position = position_jitter(width = 0.08, seed = 42)) +
  scale_fill_manual(values = c("Control" = "#0279EE", "ONFH" = "#FF9400")) +
  labs(title = "Fibrocartilage program, donor-level UCell",
       subtitle = sprintf("Donors as biological units (n=3 vs 3); Welch t-test p=%.4f", tt_u$p.value),
       x = NULL, y = "Donor mean UCell score") +
  theme_pub
save_fig(p_uc, "Fig2_module_scores_donor_UCell", 6, 5)

# ---- 4. Audit table ----
audit <- data.frame(
  figure = c("Fig2_module_scores.svg", "Fig2b_umap_feature_scores.svg",
             "R2_Fig_donor_scores.svg", "R2_Fig_donor_scores_bar.svg",
             "R2_Fig_donor_scores_corrected.svg", "Fig2_module_scores_donor_AddModuleScore.svg",
             "Fig2_module_scores_donor_UCell.svg"),
  algorithm = c("UCell", "UCell", "AddModuleScore (MS_1-5)", "AddModuleScore (MS_1-5)",
                "AddModuleScore (MS_1-5)", "AddModuleScore (MS_2)", "UCell"),
  level = c("cell (n=57,310)", "cell (n=57,310)", "donor (n=6)", "donor (n=6)",
            "donor (n=6)", "donor (n=6)", "donor (n=6)"),
  label_was_correct = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
  consistent_with_Results_3_2_text = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE),
  note = c("title states UCell; keep only if text quotes UCell statistics",
           "no algorithm label; 0-0.6 non-negative scale = UCell",
           "data are AddModuleScore but y axis said 'UCell score' - label bug",
           "same label bug", "corrected label; data unchanged",
           "matches text values 0.556 vs 1.830 (MS_2 donor means)",
           "UCell counterpart: 0.407 vs 0.706, Welch p=0.0094, d=7.43"))
write.csv(audit, file.path(TAB, "R2_figure2_scoring_audit.csv"), row.names = FALSE)
print(audit)
cat("\nDONE figure2 scoring audit\n")
