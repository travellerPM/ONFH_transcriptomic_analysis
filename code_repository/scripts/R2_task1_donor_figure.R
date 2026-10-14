#!/usr/bin/env R
# R2 Task 1: Donor-level figure
# Replace cell-level boxplots with donor-level mean scores
# Panel A: Fibrocartilage donor mean scores (boxplot + individual points, n=3/group)
# Panel B: All 5 programs donor-level mean scores (dot plot, facet)

suppressMessages({
  library(ggplot2)
  library(svglite)
  library(patchwork)
})

cat("=== R2 Task 1: Donor-Level Figure ===\n")

# --- Load donor module scores ---
dm <- read.csv("/mnt/results/phase1/tables/donor_module_scores.csv")
cat("Donor module scores:", nrow(dm), "rows\n")
cat("Columns:", paste(colnames(dm), collapse = ", "), "\n")
cat("Donors:", paste(unique(dm$donor), collapse = ", "), "\n")
cat("Programs:", paste(unique(dm$program), collapse = ", "), "\n")

# Set factor order
dm$condition <- factor(dm$condition, levels = c("Control", "ONFH"))
dm$donor <- factor(dm$donor, levels = c("Ctrl-1", "Ctrl-2", "Ctrl-3", "ONFH-1", "ONFH-2", "ONFH-3"))

theme_pub <- theme_bw() +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 9))

# --- Panel A: Fibrocartilage donor mean scores ---
cat("\n--- Panel A: Fibrocartilage ---\n")
fc <- dm[dm$program == "Fibrocartilage", ]
cat("Fibrocartilage donor scores:\n")
print(fc[, c("donor", "condition", "mean_score", "n_cells")])

# Group statistics
fc_stats <- aggregate(mean_score ~ condition, data = fc,
                       FUN = function(x) c(mean = mean(x), sd = sd(x)))
fc_stats <- cbind(fc_stats[, 1], as.data.frame(fc_stats[, 2]))
colnames(fc_stats) <- c("condition", "mean", "sd")
cat("\nGroup statistics:\n")
print(fc_stats)

# Welch t-test (n=3 per group)
tt <- t.test(mean_score ~ condition, data = fc)
cat(sprintf("\nWelch t-test: t=%.3f, df=%.2f, p=%.4f\n", tt$statistic, tt$parameter, tt$p.value))

# Effect size (Cohen's d)
d_cohen <- (fc_stats$mean[fc_stats$condition == "ONFH"] - fc_stats$mean[fc_stats$condition == "Control"]) /
            sqrt((fc_stats$sd[1]^2 + fc_stats$sd[2]^2) / 2)
cat(sprintf("Cohen's d: %.2f\n", d_cohen))

# Save stats
fc_stats_out <- data.frame(
  program = "Fibrocartilage",
  condition = fc_stats$condition,
  n_donors = 3,
  mean_score = round(fc_stats$mean, 4),
  sd_score = round(fc_stats$sd, 4),
  t_stat = round(tt$statistic, 3),
  df = round(tt$parameter, 2),
  p_value = format.pval(tt$p.value, digits = 3),
  cohens_d = round(d_cohen, 2)
)
write.csv(fc_stats_out, "/mnt/results/phase1/tables/R2_donor_scores_group_stats.csv", row.names = FALSE)

# Panel A plot
p_a <- ggplot(fc, aes(x = condition, y = mean_score, fill = condition)) +
  geom_boxplot(width = 0.4, alpha = 0.3, outlier.shape = NA) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.8,
             position = position_jitter(width = 0.1, seed = 42)) +
  geom_line(aes(group = NA), data = NULL) +  # placeholder
  scale_fill_manual(values = c("Control" = "#0279EE", "ONFH" = "#FF9400")) +
  labs(title = "Fibrocartilage remodeling program",
       subtitle = sprintf("Donor-level mean UCell score (n=3 per group); Welch t-test p=%.4f, Cohen's d=%.2f",
                           tt$p.value, d_cohen),
       x = NULL, y = "Donor mean UCell score") +
  theme_pub + theme(legend.position = "none")

# --- Panel B: All 5 programs ---
cat("\n--- Panel B: All programs ---\n")

# Compute stats for all programs
all_stats <- data.frame()
for (prog in unique(dm$program)) {
  sub <- dm[dm$program == prog, ]
  tt_p <- t.test(mean_score ~ condition, data = sub)
  ctrl_mean <- mean(sub$mean_score[sub$condition == "Control"])
  ctrl_sd <- sd(sub$mean_score[sub$condition == "Control"])
  onfh_mean <- mean(sub$mean_score[sub$condition == "ONFH"])
  onfh_sd <- sd(sub$mean_score[sub$condition == "ONFH"])
  all_stats <- rbind(all_stats, data.frame(
    program = prog,
    control_mean = ctrl_mean,
    control_sd = ctrl_sd,
    onfh_mean = onfh_mean,
    onfh_sd = onfh_sd,
    p_value = tt_p$p.value
  ))
}
cat("\nAll program stats:\n")
print(all_stats)

# Save all stats
write.csv(all_stats, "/mnt/results/phase1/tables/R2_donor_scores_all_programs_stats.csv", row.names = FALSE)

# Panel B plot
p_b <- ggplot(dm, aes(x = condition, y = mean_score, fill = condition)) +
  geom_boxplot(width = 0.4, alpha = 0.3, outlier.shape = NA) +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.7,
             position = position_jitter(width = 0.1, seed = 42)) +
  facet_wrap(~program, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Control" = "#0279EE", "ONFH" = "#FF9400")) +
  labs(title = "All remodeling programs",
       x = NULL, y = "Donor mean UCell score") +
  theme_pub + theme(legend.position = "none")

# Combine
p_combined <- p_a / p_b + plot_layout(heights = c(1, 1))

svglite::svglite("/mnt/results/phase1/figures/R2_Fig_donor_scores.svg", width = 14, height = 8)
print(p_combined)
dev.off()
ggsave("/mnt/results/phase1/figures/R2_Fig_donor_scores.png", p_combined, width = 14, height = 8, dpi = 300, bg = "white")
cat("Saved R2_Fig_donor_scores\n")

# Also save a simple bar+point version (alternative visualization)
p_bar <- ggplot(dm, aes(x = program, y = mean_score, fill = condition)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5, alpha = 0.6) +
  geom_point(aes(group = condition), position = position_dodge(width = 0.6),
             size = 3, shape = 21, color = "black", stroke = 0.7) +
  scale_fill_manual(values = c("Control" = "#0279EE", "ONFH" = "#FF9400")) +
  labs(x = NULL, y = "Donor mean UCell score",
       title = "Donor-level program scores (n=3 per group)") +
  theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
                     legend.position = "bottom", legend.title = element_blank())

svglite::svglite("/mnt/results/phase1/figures/R2_Fig_donor_scores_bar.svg", width = 10, height = 6)
print(p_bar)
dev.off()
ggsave("/mnt/results/phase1/figures/R2_Fig_donor_scores_bar.png", p_bar, width = 10, height = 6, dpi = 300, bg = "white")
cat("Saved R2_Fig_donor_scores_bar\n")

cat("\n=== R2 Task 1 Complete ===\n")
