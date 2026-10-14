#!/usr/bin/env R
# R2 Task 3: scDblFinder doublet detection on GSE311084 + GSE290411
# Reports per-sample doublet rate + before/after QC comparison

suppressMessages({
  library(Seurat)
  library(scDblFinder)
  library(Matrix)
})

cat("=== R2 Task 3: scDblFinder Doublet Detection ===\n")

# ---------------------------------------------------------------
# GSE311084
# ---------------------------------------------------------------
cat("\n--- GSE311084 ---\n")
seu <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored.rds")
seu <- JoinLayers(seu)
cat("Cells:", ncol(seu), " Genes:", nrow(seu), "\n")

# scDblFinder needs a SingleCellExperiment or counts matrix
counts_311 <- GetAssayData(seu, assay = "RNA", layer = "counts")
cat("Counts matrix:", nrow(counts_311), "x", ncol(counts_311), "\n")

# Run scDblFinder per sample (batch = donor_id)
cat("Running scDblFinder (batch = donor_id)...\n")
sce_311 <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts_311)
)
colData(sce_311)$donor_id <- seu$donor_id
colData(sce_311)$condition <- seu$condition

set.seed(42)
sce_311 <- scDblFinder(sce_311, samples = colData(sce_311)$donor_id, verbose = TRUE)

# Extract results
dbl_311 <- as.data.frame(colData(sce_311))
cat("\nDoublet results columns:", paste(grep("dbl|scDbl", colnames(dbl_311), value = TRUE), collapse = ", "), "\n")

# Per-sample summary
dbl_summary_311 <- aggregate(
  list(n_cells = dbl_311$scDblFinder.class,
       n_doublet = as.numeric(dbl_311$scDblFinder.class == "doublet")),
  by = list(donor_id = dbl_311$donor_id, condition = dbl_311$condition),
  FUN = function(x) if (is.numeric(x)) sum(x) else length(x)
)
# Fix: recount properly
dbl_summary_311 <- data.frame(
  donor_id = sort(unique(dbl_311$donor_id)),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(dbl_summary_311))) {
  d <- dbl_summary_311$donor_id[i]
  sub <- dbl_311[dbl_311$donor_id == d, ]
  dbl_summary_311$condition[i] <- unique(sub$condition)[1]
  dbl_summary_311$n_cells_before[i] <- nrow(sub)
  dbl_summary_311$n_doublet[i] <- sum(sub$scDblFinder.class == "doublet")
  dbl_summary_311$n_singlet[i] <- sum(sub$scDblFinder.class == "singlet")
  dbl_summary_311$doublet_rate[i] <- round(100 * sum(sub$scDblFinder.class == "doublet") / nrow(sub), 2)
  dbl_summary_311$mean_doublet_score[i] <- round(mean(sub$scDblFinder.score, na.rm = TRUE), 4)
}

cat("\nGSE311084 doublet rates per sample:\n")
print(dbl_summary_311)

# Overall
cat("\nOverall GSE311084 doublet rate:",
    round(100 * sum(dbl_311$scDblFinder.class == "doublet") / nrow(dbl_311), 2), "%\n")

# Save per-cell doublet calls
write.csv(dbl_311[, c("donor_id", "condition", "scDblFinder.class", "scDblFinder.score", "scDblFinder.weighted")],
          "/mnt/results/phase1/tables/R2_doublet_calls_GSE311084_per_cell.csv", row.names = TRUE)

# Save per-sample summary
write.csv(dbl_summary_311, "/mnt/results/phase1/tables/R2_doublet_rates_GSE311084.csv", row.names = FALSE)

# Before/after QC comparison
qc_before_after_311 <- data.frame(
  dataset = "GSE311084",
  donor_id = dbl_summary_311$donor_id,
  condition = dbl_summary_311$condition,
  n_cells_before_doublet_removal = dbl_summary_311$n_cells_before,
  n_doublets_removed = dbl_summary_311$n_doublet,
  n_cells_after_doublet_removal = dbl_summary_311$n_singlet,
  doublet_rate_pct = dbl_summary_311$doublet_rate
)
write.csv(qc_before_after_311, "/mnt/results/phase1/tables/R2_doublet_QC_GSE311084.csv", row.names = FALSE)

# Save the Seurat object with doublet calls added
seu$scDblFinder_class <- dbl_311$scDblFinder.class
seu$scDblFinder_score <- dbl_311$scDblFinder.score
saveRDS(seu, "/mnt/shared-workspace/shared/phase1/GSE311084_scored_doublet.rds")
cat("Saved GSE311084_scored_doublet.rds\n")

rm(seu, counts_311, sce_311, dbl_311); gc()

# ---------------------------------------------------------------
# GSE290411
# ---------------------------------------------------------------
cat("\n--- GSE290411 ---\n")
seu2 <- readRDS("/mnt/shared-workspace/shared/phase1/GSE290411_scored.rds")
seu2 <- JoinLayers(seu2)
cat("Cells:", ncol(seu2), " Genes:", nrow(seu2), "\n")

# Check per-donor cell counts — exclude donors with too few cells for scDblFinder
donor_counts <- table(seu2$donor_id)
cat("Per-donor cell counts:\n")
print(donor_counts)
# Keep only donors with >= 10 cells (scDblFinder minimum)
keep_donors <- names(donor_counts[donor_counts >= 10])
cat("Donors with >= 10 cells:", paste(keep_donors, collapse = ", "), "\n")
seu2_sub <- subset(seu2, donor_id %in% keep_donors)
cat("Cells after excluding small donors:", ncol(seu2_sub), "\n")

counts_290 <- GetAssayData(seu2_sub, assay = "RNA", layer = "counts")
sce_290 <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts_290)
)
colData(sce_290)$donor_id <- seu2_sub$donor_id
colData(sce_290)$condition <- seu2_sub$condition

set.seed(42)
sce_290 <- scDblFinder(sce_290, samples = colData(sce_290)$donor_id, verbose = TRUE)

dbl_290 <- as.data.frame(colData(sce_290))

dbl_summary_290 <- data.frame(
  donor_id = sort(unique(dbl_290$donor_id)),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(dbl_summary_290))) {
  d <- dbl_summary_290$donor_id[i]
  sub <- dbl_290[dbl_290$donor_id == d, ]
  dbl_summary_290$condition[i] <- unique(sub$condition)[1]
  dbl_summary_290$n_cells_before[i] <- nrow(sub)
  dbl_summary_290$n_doublet[i] <- sum(sub$scDblFinder.class == "doublet")
  dbl_summary_290$n_singlet[i] <- sum(sub$scDblFinder.class == "singlet")
  dbl_summary_290$doublet_rate[i] <- round(100 * sum(sub$scDblFinder.class == "doublet") / nrow(sub), 2)
  dbl_summary_290$mean_doublet_score[i] <- round(mean(sub$scDblFinder.score, na.rm = TRUE), 4)
}

cat("\nGSE290411 doublet rates per sample:\n")
print(dbl_summary_290)

cat("\nOverall GSE290411 doublet rate:",
    round(100 * sum(dbl_290$scDblFinder.class == "doublet") / nrow(dbl_290), 2), "%\n")

write.csv(dbl_290[, c("donor_id", "condition", "scDblFinder.class", "scDblFinder.score", "scDblFinder.weighted")],
          "/mnt/results/phase1/tables/R2_doublet_calls_GSE290411_per_cell.csv", row.names = TRUE)
write.csv(dbl_summary_290, "/mnt/results/phase1/tables/R2_doublet_rates_GSE290411.csv", row.names = FALSE)

qc_before_after_290 <- data.frame(
  dataset = "GSE290411",
  donor_id = dbl_summary_290$donor_id,
  condition = dbl_summary_290$condition,
  n_cells_before_doublet_removal = dbl_summary_290$n_cells_before,
  n_doublets_removed = dbl_summary_290$n_doublet,
  n_cells_after_doublet_removal = dbl_summary_290$n_singlet,
  doublet_rate_pct = dbl_summary_290$doublet_rate
)
write.csv(qc_before_after_290, "/mnt/results/phase1/tables/R2_doublet_QC_GSE290411.csv", row.names = FALSE)

seu2$scDblFinder_class <- dbl_290$scDblFinder.class
seu2$scDblFinder_score <- dbl_290$scDblFinder.score
saveRDS(seu2, "/mnt/shared-workspace/shared/phase1/GSE290411_scored_doublet.rds")
cat("Saved GSE290411_scored_doublet.rds\n")

# ---------------------------------------------------------------
# Combined summary
# ---------------------------------------------------------------
cat("\n--- Combined Doublet Summary ---\n")
combined <- rbind(
  data.frame(dataset = "GSE311084", dbl_summary_311),
  data.frame(dataset = "GSE290411", dbl_summary_290)
)
write.csv(combined, "/mnt/results/phase1/tables/R2_doublet_rates_combined.csv", row.names = FALSE)
print(combined)

cat("\n=== R2 Task 3 Complete ===\n")
