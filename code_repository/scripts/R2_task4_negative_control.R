#!/usr/bin/env R
# R2 Task 4: Negative control validation
# 1000 random gene sets matched for gene number + expression level
# Compare observed fibrocartilage score to empirical null
# Also compare with Hallmark EMT + generic ECM/fibroblast signatures

suppressMessages({
  library(Seurat)
  library(UCell)
  library(msigdbr)
})

cat("=== R2 Task 4: Negative Control Validation ===\n")

# --- Load scored object ---
seu <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored.rds")
seu <- JoinLayers(seu)
cat("Cells:", ncol(seu), " Genes:", nrow(seu), "\n")

# --- Load frozen signatures ---
frozen <- read.csv("/mnt/results/phase1/tables/frozen_signatures.csv")
program_names <- c("Homeostatic", "Fibrocartilage", "Inflammatory",
                   "Pre_hypertrophic", "Hypertrophic_osteogenic")
sig_lists <- lapply(program_names, function(p) {
  frozen$gene[frozen$signature_name == p & frozen$gene %in% rownames(seu)]
})
names(sig_lists) <- program_names
cat("Signature sizes (present in dataset):\n")
for (nm in program_names) cat("  ", nm, ":", length(sig_lists[[nm]]), "\n")

# --- Compute mean expression for expression matching ---
cat("\nComputing mean expression for expression-bin matching...\n")
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
# Use log1p of mean counts as expression level proxy
gene_means <- Matrix::rowMeans(counts)
gene_expr_log <- log1p(gene_means)
cat("Genes with mean > 0:", sum(gene_means > 0), "\n")

# Create expression bins (deciles among expressed genes)
expressed_genes <- names(gene_expr_log)[gene_means > 0]
expr_bins <- cut(gene_expr_log[expressed_genes], breaks = 10, labels = FALSE)
names(expr_bins) <- expressed_genes
cat("Expression bins:\n")
print(table(expr_bins))

# --- Generate 1000 random gene sets per program (expression-matched) ---
cat("\nGenerating 1000 expression-matched random gene sets per program...\n")
set.seed(42)
N_RANDOM <- 1000

# For each program, find which bins its genes fall into
random_sets <- list()
for (prog in program_names) {
  genes <- sig_lists[[prog]]
  # Find bins of the program genes
  prog_bins <- expr_bins[genes]
  prog_bins <- prog_bins[!is.na(prog_bins)]
  bin_table <- table(prog_bins)

  # Sample from the same bins
  sets <- list()
  for (i in 1:N_RANDOM) {
    sampled <- unlist(lapply(names(bin_table), function(b) {
      pool <- names(expr_bins)[expr_bins == as.numeric(b)]
      sample(pool, size = bin_table[b], replace = FALSE)
    }))
    sets[[i]] <- unique(sampled)
  }
  random_sets[[prog]] <- sets
  cat("  ", prog, ":", length(genes), "genes, bins:", paste(names(bin_table), collapse = ","), "\n")
}

# --- Score all random sets + observed signatures with UCell ---
cat("\nScoring with UCell (this may take a while)...\n")

# Build a single feature list: observed + all random sets
all_features <- list()
for (prog in program_names) {
  all_features[[paste0("OBS_", prog)]] <- sig_lists[[prog]]
  for (i in 1:N_RANDOM) {
    all_features[[paste0("RND_", prog, "_", i)]] <- random_sets[[prog]][[i]]
  }
}
cat("Total gene sets to score:", length(all_features), "\n")

# Run UCell on a stratified subset for speed (10k cells)
set.seed(42)
meta <- seu@meta.data
cells_sub <- unlist(lapply(split(colnames(seu), meta$condition), function(x) sample(x, min(5000, length(x)))))
seu_sub <- subset(seu, cells = cells_sub)
cat("Scoring on", ncol(seu_sub), "cells (stratified subset)\n")

seu_scored <- AddModuleScore_UCell(seu_sub, features = all_features, name = "NC_")
cat("UCell scoring complete\n")

# --- Extract results ---
meta_nc <- seu_scored@meta.data
# UCell appends name as suffix: {feature_name}NC_
nc_cols <- grep("NC_$", colnames(meta_nc), value = TRUE)
cat("Score columns:", length(nc_cols), "\n")

# --- Compute null distribution and empirical p-values ---
cat("\n--- Negative Control Results ---\n")
nc_results <- data.frame()
null_distributions <- list()

for (prog in program_names) {
  obs_col <- paste0("OBS_", prog, "NC_")
  obs_scores <- meta_nc[[obs_col]]

  # Get random scores for this program
  rnd_cols <- grep(paste0("^RND_", prog, "_.*NC_$"), colnames(meta_nc), value = TRUE)
  rnd_mean_scores <- colMeans(meta_nc[, rnd_cols, drop = FALSE])

  # Null distribution: mean score across cells for each random set
  null_mean <- mean(rnd_mean_scores)
  null_sd <- sd(rnd_mean_scores)

  # Observed: mean score across cells
  obs_mean <- mean(obs_scores)

  # Empirical p-value (two-sided: |observed - null| relative to null distribution)
  n_greater <- sum(rnd_mean_scores >= obs_mean)
  n_lesser <- sum(rnd_mean_scores <= obs_mean)
  emp_p_upper <- (1 + n_greater) / (1 + N_RANDOM)
  emp_p_lower <- (1 + n_lesser) / (1 + N_RANDOM)
  emp_p <- min(emp_p_upper, emp_p_lower) * 2  # two-sided
  emp_p <- min(emp_p, 1)  # cap at 1

  # Z-score relative to null
  z_score <- (obs_mean - null_mean) / null_sd

  # Disease effect: ONFH vs Control difference
  obs_diff <- mean(obs_scores[meta_nc$condition == "ONFH"]) -
              mean(obs_scores[meta_nc$condition == "Control"])

  # Random sets disease effect
  rnd_diffs <- sapply(rnd_cols, function(rc) {
    mean(meta_nc[[rc]][meta_nc$condition == "ONFH"]) -
    mean(meta_nc[[rc]][meta_nc$condition == "Control"])
  })
  null_diff_mean <- mean(rnd_diffs)
  null_diff_sd <- sd(rnd_diffs)
  n_diff_greater <- sum(rnd_diffs >= obs_diff)
  n_diff_lesser <- sum(rnd_diffs <= obs_diff)
  emp_p_diff_upper <- (1 + n_diff_greater) / (1 + N_RANDOM)
  emp_p_diff_lower <- (1 + n_diff_lesser) / (1 + N_RANDOM)
  emp_p_diff <- min(min(emp_p_diff_upper, emp_p_diff_lower) * 2, 1)  # two-sided
  z_diff <- (obs_diff - null_diff_mean) / null_diff_sd

  nc_results <- rbind(nc_results, data.frame(
    program = prog,
    n_genes = length(sig_lists[[prog]]),
    observed_mean_score = round(obs_mean, 4),
    null_mean_score = round(null_mean, 4),
    null_sd = round(null_sd, 4),
    z_score = round(z_score, 2),
    empirical_p_score = emp_p,
    observed_disease_diff = round(obs_diff, 4),
    null_mean_disease_diff = round(null_diff_mean, 4),
    null_sd_disease_diff = round(null_diff_sd, 4),
    z_score_disease_diff = round(z_diff, 2),
    empirical_p_disease_diff = emp_p_diff
  ))

  # Store null distribution for plotting
  null_distributions[[prog]] <- rnd_mean_scores

  cat(sprintf("  %s: obs=%.4f, null=%.4f±%.4f, z=%.2f, p=%.4f | disease_diff=%.4f, null_diff=%.4f±%.4f, z=%.2f, p=%.4f\n",
      prog, obs_mean, null_mean, null_sd, z_score, emp_p,
      obs_diff, null_diff_mean, null_diff_sd, z_diff, emp_p_diff))
}

# Save results
write.csv(nc_results, "/mnt/results/phase1/tables/R2_negative_control_results.csv", row.names = FALSE)

# Save null distributions
null_df <- do.call(cbind, lapply(names(null_distributions), function(p) {
  data.frame(program = p, random_set_id = 1:N_RANDOM, mean_score = null_distributions[[p]])
}))
write.csv(null_df, "/mnt/results/phase1/tables/R2_negative_control_null_distribution.csv", row.names = FALSE)

# --- Comparison with Hallmark EMT + generic ECM/fibroblast signatures ---
cat("\n--- Comparison with EMT / ECM / Fibroblast signatures ---\n")

# Get Hallmark EMT from msigdbr
msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
emt_genes <- msig_h$gene_symbol[msig_h$gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
emt_genes <- emt_genes[emt_genes %in% rownames(seu)]
cat("Hallmark EMT genes present:", length(emt_genes), "\n")

# Generic ECM/fibroblast signature (from MSigDB C8 or custom)
# Use GO:COLLAGEN_FIBRIL_ORGANIZATION + GO:EXTRACELLULAR_MATRIX_ORGANIZATION as generic ECM
msig_c5 <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
ecm_genes1 <- msig_c5$gene_symbol[msig_c5$gs_name == "GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION"]
ecm_genes2 <- msig_c5$gene_symbol[msig_c5$gs_name == "GOBP_COLLAGEN_FIBRIL_ORGANIZATION"]
ecm_genes <- unique(c(ecm_genes1, ecm_genes2))
ecm_genes <- ecm_genes[ecm_genes %in% rownames(seu)]
cat("Generic ECM genes present:", length(ecm_genes), "\n")

# Generic fibroblast signature (from MSigDB C8: FIBROBLAST)
fib_genes <- msig_c5$gene_symbol[grepl("FIBROBLAST_PROLIFERATION", msig_c5$gs_name)]
fib_genes <- unique(fib_genes)[unique(fib_genes) %in% rownames(seu)]
if (length(fib_genes) < 5) {
  # Fallback: use a simple collagen/fibroblast marker set
  fib_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3",
                  "DCN", "LUM", "POSTN", "THBS2", "TIMP1", "BGN")
  fib_genes <- fib_genes[fib_genes %in% rownames(seu)]
}
cat("Generic fibroblast genes present:", length(fib_genes), "\n")

# Score these comparison signatures
comparison_sets <- list(
  Hallmark_EMT = emt_genes,
  Generic_ECM = ecm_genes,
  Generic_Fibroblast = fib_genes,
  Fibrocartilage_program = sig_lists[["Fibrocartilage"]]
)

seu_comp <- AddModuleScore_UCell(seu_sub, features = comparison_sets, name = "COMP_")
meta_comp <- seu_comp@meta.data

comp_results <- data.frame()
for (nm in names(comparison_sets)) {
  col <- paste0(nm, "COMP_")
  scores <- meta_comp[[col]]
  obs_diff <- mean(scores[meta_comp$condition == "ONFH"]) -
              mean(scores[meta_comp$condition == "Control"])

  # Wilcoxon rank-sum test (donor-level would be better but n=3)
  wt <- wilcox.test(scores[meta_comp$condition == "ONFH"],
                    scores[meta_comp$condition == "Control"])

  comp_results <- rbind(comp_results, data.frame(
    signature = nm,
    n_genes = length(comparison_sets[[nm]]),
    control_mean = round(mean(scores[meta_comp$condition == "Control"]), 4),
    onfh_mean = round(mean(scores[meta_comp$condition == "ONFH"]), 4),
    disease_diff = round(obs_diff, 4),
    wilcoxon_p = format.pval(wt$p.value, digits = 3)
  ))
}
cat("\nComparison signatures:\n")
print(comp_results)

write.csv(comp_results, "/mnt/results/phase1/tables/R2_comparison_signatures.csv", row.names = FALSE)

# --- Figure: Null distribution + observed position ---
cat("\n--- Creating Negative Control Figure ---\n")
library(ggplot2)
library(svglite)
library(patchwork)

theme_pub <- theme_bw() +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 9))

# Build plot data
plot_data <- do.call(rbind, lapply(names(null_distributions), function(p) {
  data.frame(program = p, mean_score = null_distributions[[p]], type = "Null (n=1000)")
}))
obs_data <- nc_results[, c("program", "observed_mean_score")]
colnames(obs_data) <- c("program", "mean_score")
obs_data$type <- "Observed"

# Panel A: Null distributions with observed position
p_a <- ggplot(plot_data, aes(x = mean_score, fill = program)) +
  geom_density(alpha = 0.4, color = "grey40", linewidth = 0.3) +
  geom_vline(data = obs_data, aes(xintercept = mean_score, color = program),
             linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~program, scales = "free") +
  scale_fill_manual(values = rep("grey80", 5), guide = "none") +
  scale_color_manual(values = c("Homeostatic" = "#75A025", "Fibrocartilage" = "#FF9400",
                                 "Inflammatory" = "#FD9BED", "Pre_hypertrophic" = "#0279EE",
                                 "Hypertrophic_osteogenic" = "#E9ED4C"), guide = "none") +
  labs(x = "Mean UCell score (across cells)", y = "Density",
       title = "Negative control: expression-matched random gene sets (n=1000)") +
  theme_pub

# Panel B: Disease effect (ONFH - Control) for observed vs null
diff_data <- do.call(rbind, lapply(program_names, function(p) {
  rnd_cols <- grep(paste0("^RND_", p, "_.*NC_$"), colnames(meta_nc), value = TRUE)
  rnd_diffs <- sapply(rnd_cols, function(rc) {
    mean(meta_nc[[rc]][meta_nc$condition == "ONFH"]) -
    mean(meta_nc[[rc]][meta_nc$condition == "Control"])
  })
  obs_diff <- nc_results$observed_disease_diff[nc_results$program == p]
  data.frame(program = p, diff = c(rnd_diffs, obs_diff),
             type = c(rep("Null", N_RANDOM), "Observed"))
}))

p_b <- ggplot(diff_data, aes(x = program, y = diff, fill = type)) +
  geom_boxplot(data = subset(diff_data, type == "Null"), alpha = 0.5, outlier.size = 0.5) +
  geom_point(data = subset(diff_data, type == "Observed"),
             shape = 21, size = 4, fill = "#FF9400", color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("Null" = "grey80", "Observed" = "#FF9400")) +
  coord_flip() +
  labs(x = NULL, y = "Disease effect (ONFH - Control UCell score)",
       title = "Observed vs null disease effect") +
  theme_pub + theme(legend.position = "bottom", legend.title = element_blank())

p_nc <- p_a / p_b + plot_layout(heights = c(1, 1))

svglite::svglite("/mnt/results/phase1/figures/R2_Fig_negative_control.svg", width = 12, height = 10)
print(p_nc)
dev.off()
ggsave("/mnt/results/phase1/figures/R2_Fig_negative_control.png", p_nc, width = 12, height = 10, dpi = 300, bg = "white")
cat("Saved R2_Fig_negative_control\n")

cat("\n=== R2 Task 4 Complete ===\n")
