#!/usr/bin/env R
# R1 #2: Random-effects meta-analysis (DerSimonian-Laird)
# Combines GSE74089 + GSE345158 ssGSEA signature scores
# Reports combined effect, CI, p, Q, I², tau² + forest plot

suppressMessages({
  library(metafor)
  library(ggplot2)
  library(svglite)
})

cat("=== R1 #2: Random-Effects Meta-Analysis ===\n")

# --- Load validation signature scores ---
v74089 <- read.csv("/mnt/results/phase1/tables/validation_GSE74089_signature_scores.csv")
v345158 <- read.csv("/mnt/results/phase1/tables/validation_GSE345158_signature_scores.csv")

cat("GSE74089 signatures:", nrow(v74089), "\n")
cat("GSE345158 signatures:", nrow(v345158), "\n")
cat("\nGSE74089 columns:", paste(colnames(v74089), collapse = ", "), "\n")
cat("GSE345158 columns:", paste(colnames(v345158), collapse = ", "), "\n")

# --- Compute effect sizes and standard errors ---
# Effect size: ssGSEA diff (disease - control)
# SE: diff / t_stat (from the t-test)
# For GSE74089: n_NFH=4, n_Ctrl=4, df=6
# For GSE345158: n_SONFH=5, n_Ctrl=4, df=7

compute_se <- function(diff, t_stat, n1, n2) {
  # For a two-sample t-test: SE = diff / t
  se <- abs(diff / t_stat)
  return(se)
}

# GSE74089: 4 vs 4
v74089$dataset <- "GSE74089"
v74089$se <- compute_se(v74089$diff, v74089$t_stat, 4, 4)
v74089$n_disease <- 4
v74089$n_control <- 4

# GSE345158: 5 vs 4
v345158$dataset <- "GSE345158"
v345158$se <- compute_se(v345158$diff, v345158$t_stat, 5, 4)
v345158$n_disease <- 5
v345158$n_control <- 4

# Combine
combined <- rbind(
  v74089[, c("signature", "dataset", "diff", "se", "p_value", "n_disease", "n_control")],
  v345158[, c("signature", "dataset", "diff", "se", "p_value", "n_disease", "n_control")]
)
colnames(combined)[colnames(combined) == "diff"] <- "effect_size"
colnames(combined)[colnames(combined) == "p_value"] <- "p_value"

cat("\nCombined effect sizes:\n")
print(combined)

# --- Run random-effects meta-analysis for each signature ---
signatures <- unique(combined$signature)
meta_results <- data.frame()

for (sig in signatures) {
  sub <- combined[combined$signature == sig, ]
  if (nrow(sub) < 2) next

  # Random-effects meta-analysis (DerSimonian-Laird)
  rma_fit <- rma(yi = effect_size, sei = se, data = sub,
                  method = "DL", test = "z")

  meta_results <- rbind(meta_results, data.frame(
    signature = sig,
    n_datasets = nrow(sub),
    combined_effect = round(rma_fit$beta[1], 4),
    se_combined = round(rma_fit$se, 4),
    ci_lower = round(rma_fit$ci.lb, 4),
    ci_upper = round(rma_fit$ci.ub, 4),
    z_value = round(rma_fit$zval, 3),
    p_value = format.pval(rma_fit$pval, digits = 3),
    Q_statistic = round(rma_fit$QE, 3),
    Q_df = rma_fit$k - 1,
    Q_pvalue = format.pval(rma_fit$QEp, digits = 3),
    I2 = round(rma_fit$I2, 1),
    tau2 = round(rma_fit$tau2, 6),
    method = "DerSimonian-Laird random-effects"
  ))
}

cat("\n--- Meta-Analysis Results ---\n")
print(meta_results)

write.csv(meta_results, "/mnt/results/phase1/tables/R1_meta_analysis_results.csv", row.names = FALSE)

# --- Forest plot ---
cat("\n--- Creating Forest Plot ---\n")

# Prepare data for forest plot
forest_data <- combined
forest_data$dataset_label <- ifelse(forest_data$dataset == "GSE74089",
                                     "GSE74089\n(hip cartilage, microarray)",
                                     "GSE345158\n(femoral head, RNA-seq)")

# Add combined estimate
for (i in seq_len(nrow(meta_results))) {
  sig <- meta_results$signature[i]
  forest_data <- rbind(forest_data, data.frame(
    signature = sig,
    dataset = "Combined",
    effect_size = meta_results$combined_effect[i],
    se = meta_results$se_combined[i],
    p_value = as.numeric(meta_results$p_value[i]),
    n_disease = NA,
    n_control = NA,
    dataset_label = "Combined (RE)"
  ))
}

# Compute CI for each row
forest_data$ci_lower <- forest_data$effect_size - 1.96 * forest_data$se
forest_data$ci_upper <- forest_data$effect_size + 1.96 * forest_data$se

# Order signatures by combined effect
sig_order <- meta_results$signature[order(meta_results$combined_effect)]
forest_data$signature <- factor(forest_data$signature, levels = rev(sig_order))
forest_data$dataset <- factor(forest_data$dataset,
                               levels = c("GSE74089", "GSE345158", "Combined"))

# Forest plot
p_forest <- ggplot(forest_data, aes(x = effect_size, y = dataset, color = dataset)) +
  geom_point(size = 3, shape = 18) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.15, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  facet_wrap(~signature, ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("GSE74089" = "#0279EE", "GSE345158" = "#75A025",
                                 "Combined" = "#FF9400")) +
  labs(x = "ssGSEA effect size (disease - control)", y = NULL,
       title = "Random-effects meta-analysis: External validation signatures",
       subtitle = "DerSimonian-Laird; error bars = 95% CI") +
  theme_bw() +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_blank())

svglite::svglite("/mnt/results/phase1/figures/R1_Fig_meta_forest.svg", width = 10, height = 12)
print(p_forest)
dev.off()
ggsave("/mnt/results/phase1/figures/R1_Fig_meta_forest.png", p_forest, width = 10, height = 12, dpi = 300, bg = "white")
cat("Saved R1_Fig_meta_forest\n")

# --- Also create a simple summary table for the response ---
cat("\n--- Summary for Response ---\n")
for (i in seq_len(nrow(meta_results))) {
  cat(sprintf("  %s: combined=%.4f [%.4f, %.4f], p=%s, I²=%.1f%%, Q=%s (p=%s)\n",
      meta_results$signature[i],
      meta_results$combined_effect[i],
      meta_results$ci_lower[i],
      meta_results$ci_upper[i],
      meta_results$p_value[i],
      meta_results$I2[i],
      meta_results$Q_statistic[i],
      meta_results$Q_pvalue[i]))
}

cat("\n=== R1 #2 Meta-Analysis Complete ===\n")
