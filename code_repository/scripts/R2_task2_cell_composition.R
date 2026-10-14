#!/usr/bin/env R
# R2 Task 2: Cell composition per sample
# Table 1: before re-annotation (phase0) vs after re-annotation (lineage_call)
# Table 2: canonical chondrocyte markers (SOX9, ACAN, COL2A1) in fibrocartilage-dominant cells

suppressMessages({
  library(Seurat)
  library(jsonlite)
  library(ggplot2)
  library(svglite)
})

cat("=== R2 Task 2: Cell Composition Per Sample ===\n")

# --- 1. Load phase0 results (before re-annotation) ---
phase0 <- fromJSON("/mnt/shared-workspace/shared/onfh_data/phase0_celltype_results.json", simplifyDataFrame = FALSE)

# Build before-re-annotation table
before_list <- list()
for (i in seq_along(phase0)) {
  s <- phase0[[i]]
  ct <- s$cell_type_counts
  if (is.data.frame(ct)) {
    for (cn in colnames(ct)) {
      before_list[[length(before_list) + 1]] <- data.frame(
        dataset = s$dataset,
        sample_id = s$sample_id,
        donor_id = s$donor_id,
        disease_status = s$disease_status,
        n_cells_raw = s$n_cells_raw,
        n_cells_after_QC = s$n_cells_filtered,
        cell_type_phase0 = cn,
        n_cells_phase0 = as.integer(ct[1, cn])
      )
    }
  } else if (is.list(ct)) {
    for (nm in names(ct)) {
      before_list[[length(before_list) + 1]] <- data.frame(
        dataset = s$dataset,
        sample_id = s$sample_id,
        donor_id = s$donor_id,
        disease_status = s$disease_status,
        n_cells_raw = s$n_cells_raw,
        n_cells_after_QC = s$n_cells_filtered,
        cell_type_phase0 = nm,
        n_cells_phase0 = as.integer(ct[[nm]])
      )
    }
  }
}
before_df <- do.call(rbind, before_list)
cat("Before re-annotation:", nrow(before_df), "rows\n")
print(head(before_df, 10))

# --- 2. Load Seurat object (after re-annotation) ---
cat("\nLoading GSE311084 scored object...\n")
seu <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored.rds")
seu <- JoinLayers(seu)
meta <- seu@meta.data

# After re-annotation: lineage_call per donor
after_df <- data.frame(
  donor_id = names(table(meta$donor_id)),
  n_cells_lineage_call = as.integer(table(meta$donor_id)),
  cartilage_lineage = as.integer(table(meta$donor_id, meta$lineage_call)[, "cartilage_lineage"]),
  cartilage_fibrocartilage_candidate = as.integer(table(meta$donor_id, meta$lineage_call)[, "cartilage_fibrocartilage_candidate"])
)
cat("\nAfter re-annotation:\n")
print(after_df)

# --- 3. Merge before/after ---
# For GSE311084 only (the discovery dataset)
before_311 <- before_df[before_df$dataset == "GSE311084", ]
cat("\nGSE311084 before re-annotation:\n")
print(before_311)

# Reshape: one row per donor, columns for each phase0 cell type
before_wide <- reshape(before_311[, c("donor_id", "disease_status", "n_cells_raw", "n_cells_after_QC", "cell_type_phase0", "n_cells_phase0")],
                        idvar = c("donor_id", "disease_status", "n_cells_raw", "n_cells_after_QC"),
                        timevar = "cell_type_phase0", direction = "wide")
colnames(before_wide) <- gsub("n_cells_phase0\\.", "n_", colnames(before_wide))
# Fill NA with 0
before_wide[is.na(before_wide)] <- 0
cat("\nBefore re-annotation (wide):\n")
print(before_wide)

# Create composition table
comp_table <- data.frame(
  sample_id = before_311$sample_id[!duplicated(before_311$donor_id)],
  donor_id = before_wide$donor_id,
  condition = before_wide$disease_status,
  n_cells_raw = before_wide$n_cells_raw,
  n_cells_after_QC = before_wide$n_cells_after_QC
)

# Add phase0 cell type columns dynamically
phase0_cols <- grep("^n_", colnames(before_wide), value = TRUE)
for (pc in phase0_cols) {
  comp_table[[paste0(pc, "_phase0")]] <- before_wide[[pc]]
}

# Add after-re-annotation info
for (i in seq_len(nrow(comp_table))) {
  d <- comp_table$donor_id[i]
  idx <- which(after_df$donor_id == d)
  if (length(idx) > 0) {
    comp_table$n_cells_final[i] <- after_df$n_cells_lineage_call[idx]
    comp_table$n_cartilage_lineage[i] <- after_df$cartilage_lineage[idx]
    comp_table$n_fibrocartilage_candidate[i] <- after_df$cartilage_fibrocartilage_candidate[idx]
  }
}

cat("\n--- Cell Composition Per Sample (GSE311084) ---\n")
print(comp_table)
write.csv(comp_table, "/mnt/results/phase1/tables/R2_cell_composition_per_sample.csv", row.names = FALSE)

# --- 4. Canonical chondrocyte markers in fibrocartilage-dominant cells ---
cat("\n--- Canonical Chondrocyte Markers ---\n")
canonical_markers <- c("SOX9", "ACAN", "COL2A1")
canonical_markers <- canonical_markers[canonical_markers %in% rownames(seu)]
cat("Markers present:", paste(canonical_markers, collapse = ", "), "\n")

# Get expression (log-normalized)
expr <- GetAssayData(seu, assay = "RNA", layer = "data")
marker_expr <- as.matrix(expr[canonical_markers, ])
cat("Marker expression matrix:", nrow(marker_expr), "x", ncol(marker_expr), "\n")

# Group by lineage_call and dominant_program
lineage <- meta$lineage_call
dom_prog <- meta$dominant_program
condition <- meta$condition
donor <- meta$donor_id

# Fibrocartilage-dominant vs other
fc_dominant <- dom_prog == "Fibrocartilage"
cat("\nFibrocartilage-dominant cells:", sum(fc_dominant), "\n")
cat("Other cells:", sum(!fc_dominant), "\n")

marker_summary <- data.frame()
for (g in canonical_markers) {
  expr_g <- marker_expr[g, ]

  # By lineage_call
  for (lc in unique(lineage)) {
    idx <- lineage == lc
    marker_summary <- rbind(marker_summary, data.frame(
      marker = g, grouping = "lineage_call", group = lc,
      n_cells = sum(idx),
      mean_expr = round(mean(expr_g[idx]), 3),
      pct_expressing = round(100 * sum(expr_g[idx] > 0) / sum(idx), 1)
    ))
  }

  # By dominant_program (Fibrocartilage vs others)
  for (dp in c("Fibrocartilage", "Homeostatic", "Inflammatory", "Pre_hypertrophic", "Hypertrophic_osteogenic", "Undetermined")) {
    idx <- dom_prog == dp
    if (sum(idx) == 0) next
    marker_summary <- rbind(marker_summary, data.frame(
      marker = g, grouping = "dominant_program", group = dp,
      n_cells = sum(idx),
      mean_expr = round(mean(expr_g[idx]), 3),
      pct_expressing = round(100 * sum(expr_g[idx] > 0) / sum(idx), 1)
    ))
  }

  # By condition
  for (cond in c("Control", "ONFH")) {
    idx <- condition == cond
    marker_summary <- rbind(marker_summary, data.frame(
      marker = g, grouping = "condition", group = cond,
      n_cells = sum(idx),
      mean_expr = round(mean(expr_g[idx]), 3),
      pct_expressing = round(100 * sum(expr_g[idx] > 0) / sum(idx), 1)
    ))
  }
}

cat("\nCanonical marker summary:\n")
print(marker_summary)
write.csv(marker_summary, "/mnt/results/phase1/tables/R2_canonical_markers_summary.csv", row.names = FALSE)

# --- 5. Figure: canonical markers in fibrocartilage-dominant cells ---
cat("\n--- Creating Canonical Marker Figure ---\n")

# Build plotting data
plot_data <- data.frame(
  marker = rep(canonical_markers, each = ncol(marker_expr)),
  cell = rep(colnames(marker_expr), times = length(canonical_markers)),
  expression = as.vector(t(marker_expr)),
  lineage_call = rep(lineage, times = length(canonical_markers)),
  dominant_program = rep(dom_prog, times = length(canonical_markers)),
  condition = rep(condition, times = length(canonical_markers)),
  stringsAsFactors = FALSE
)

# Subsample for plotting (10k cells)
set.seed(42)
plot_sub <- plot_data[sample(nrow(plot_data), min(10000, nrow(plot_data))), ]

# Add fibrocartilage-dominant label
plot_sub$fc_status <- ifelse(plot_sub$dominant_program == "Fibrocartilage", "Fibrocartilage-dominant", "Other")

p_markers <- ggplot(plot_sub, aes(x = fc_status, y = expression, fill = fc_status)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, alpha = 0.5, outlier.size = 0.3) +
  facet_wrap(~marker, scales = "free_y") +
  scale_fill_manual(values = c("Fibrocartilage-dominant" = "#FF9400", "Other" = "#0279EE")) +
  labs(x = NULL, y = "Log-normalized expression",
       title = "Canonical chondrocyte markers in fibrocartilage-dominant vs other cells") +
  theme_bw() +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.title = element_blank())

svglite::svglite("/mnt/results/phase1/figures/R2_Fig_canonical_markers.svg", width = 10, height = 5)
print(p_markers)
dev.off()
ggsave("/mnt/results/phase1/figures/R2_Fig_canonical_markers.png", p_markers, width = 10, height = 5, dpi = 300, bg = "white")
cat("Saved R2_Fig_canonical_markers\n")

cat("\n=== R2 Task 2 Complete ===\n")
