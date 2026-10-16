#!/usr/bin/env R
# =============================================================================
# R2 Task 2 (revision round 2): Per-sample canonical marker analysis
# =============================================================================
# Reviewer 2: report SOX9 / ACAN / COL2A1 detection for EACH ONFH sample,
# compared with non-fibrocartilage cartilage-lineage cells.
#
# Comparator definitions (identical to R2_task3_doublet_cartilage_lineage_comparator.R):
#   A  Fibrocartilage-dominant : dominant_program == "Fibrocartilage"
#   B  Other cartilage-lineage : lineage_call == "cartilage_lineage"
#                                & dominant_program != "Fibrocartilage"
#
# Statistics: descriptive only (n, % expressing, mean/median log-normalized
# expression over ALL cells in the group, zeros included - Seurat DotPlot
# convention). No per-sample significance tests (small n, unstable).
# Empty group -> NA with reason. Seed 42. Existing results preserved.
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

cat("=== R2 Task 2 (round 2): per-sample canonical markers ===\n")
cat("R version:", R.version.string, "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("Run started:", format(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load scored object (no doublet filtering here - marker detection analysis)
# -----------------------------------------------------------------------------
seu <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored.rds")
if (is(seu[["RNA"]], "Assay5")) seu <- JoinLayers(seu)
meta <- seu@meta.data
cat("Cells:", ncol(seu), "\n")

markers <- c("SOX9", "ACAN", "COL2A1")
stopifnot(all(markers %in% rownames(seu)))
expr <- GetAssayData(seu, assay = "RNA", layer = "data")[markers, ]
cat("Expression matrix extracted:", nrow(expr), "markers x", ncol(expr), "cells\n")
cat("Data layer range:", round(range(expr), 3), "(log-normalized)\n\n")

# -----------------------------------------------------------------------------
# 2. Group definitions (same as Task 1)
# -----------------------------------------------------------------------------
meta$is_fibro_dominant <- meta$dominant_program == "Fibrocartilage"
meta$grp <- NA
meta$grp[meta$is_fibro_dominant] <- "A_Fibrocartilage_dominant"
meta$grp[!meta$is_fibro_dominant & meta$lineage_call == "cartilage_lineage"] <- "B_Other_cartilage_lineage"
meta$grp[!meta$is_fibro_dominant & meta$lineage_call == "cartilage_fibrocartilage_candidate"] <- "C_Other_fibrocartilage_candidate"

donors <- c("Ctrl-1", "Ctrl-2", "Ctrl-3", "ONFH-1", "ONFH-2", "ONFH-3")
cond_map <- setNames(c("Control", "Control", "Control", "ONFH", "ONFH", "ONFH"), donors)
sample_map <- setNames(c("A-01", "A-03", "A-07", "A-02", "A-06", "A-08"), donors)

# -----------------------------------------------------------------------------
# 3. Per-sample x marker x group detection
# -----------------------------------------------------------------------------
rows <- list()
idx <- 1
for (d in donors) {
  for (g in c("A_Fibrocartilage_dominant", "B_Other_cartilage_lineage")) {
    sel <- meta$donor_id == d & meta$grp == g
    n <- sum(sel)
    for (mk in markers) {
      if (n == 0) {
        rows[[idx]] <- data.frame(sample_id = sample_map[[d]], donor_id = d,
          condition = cond_map[[d]], marker = mk, group_definition = g,
          n_cells = 0L, n_expressing = NA_integer_, percent_expressing = NA_real_,
          mean_log_normalized_expression = NA_real_,
          median_log_normalized_expression = NA_real_,
          note = "no cells in this group for this sample")
      } else {
        cells <- rownames(meta)[sel]
        v <- as.numeric(expr[mk, cells])
        ne <- sum(v > 0)
        rows[[idx]] <- data.frame(sample_id = sample_map[[d]], donor_id = d,
          condition = cond_map[[d]], marker = mk, group_definition = g,
          n_cells = n, n_expressing = ne, percent_expressing = 100 * ne / n,
          mean_log_normalized_expression = mean(v),
          median_log_normalized_expression = median(v),
          note = if (n < 10) paste0("small group (n=", n, "); descriptive value unstable") else "")
      }
      idx <- idx + 1
    }
  }
}
det <- do.call(rbind, rows)
cat("== Per-sample detection (all samples) ==\n"); print(det, digits = 3)

onfh_det <- det[det$condition == "ONFH", ]
cat("\n== ONFH samples only ==", nrow(onfh_det), "rows\n")

write.csv(onfh_det, file.path(OUT_TAB, "R2_canonical_markers_per_ONFH_sample.csv"), row.names = FALSE)
write.csv(det, file.path(OUT_TAB, "R2_canonical_markers_per_sample_all.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Per-sample composition (main table)
# -----------------------------------------------------------------------------
comp <- do.call(rbind, lapply(donors, function(d) {
  sel <- meta$donor_id == d
  data.frame(sample_id = sample_map[[d]], donor_id = d, condition = cond_map[[d]],
             n_cells_total = sum(sel),
             n_cartilage_lineage = sum(sel & meta$lineage_call == "cartilage_lineage"),
             n_fibrocartilage_candidate = sum(sel & meta$lineage_call == "cartilage_fibrocartilage_candidate"),
             n_fibrocartilage_dominant = sum(sel & meta$is_fibro_dominant),
             pct_fibrocartilage_dominant = 100 * sum(sel & meta$is_fibro_dominant) / sum(sel))
}))
cat("\n== Per-sample composition ==\n"); print(comp, digits = 3)
write.csv(comp, file.path(OUT_TAB, "R2_cell_composition_main_table.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Figure: dot plot, marker x sample x group
#    size = % expressing, color = mean log-normalized expression
# -----------------------------------------------------------------------------
pd <- det
pd$group_definition <- factor(pd$group_definition,
                              levels = c("A_Fibrocartilage_dominant", "B_Other_cartilage_lineage"),
                              labels = c("Fibrocartilage-dominant", "Other cartilage-lineage"))
pd$donor_id <- factor(pd$donor_id, levels = donors)
pd$marker <- factor(pd$marker, levels = markers)
pd$cond_lab <- factor(pd$condition, levels = c("Control", "ONFH"))

svg_ok <- requireNamespace("svglite", quietly = TRUE)
save_fig <- function(p, name, w = 9, h = 3.8) {
  if (svg_ok) svglite::svglite(file.path(OUT_FIG, paste0(name, ".svg")), width = w, height = h)
  else svg(file.path(OUT_FIG, paste0(name, ".svg")), width = w, height = h)
  print(p); dev.off()
  ggsave(file.path(OUT_FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

p <- ggplot(pd, aes(x = donor_id, y = marker)) +
  geom_point(aes(size = percent_expressing, color = mean_log_normalized_expression,
                 fill = mean_log_normalized_expression),
             shape = 21, stroke = 0.4, alpha = 0.9) +
  geom_text(aes(label = ifelse(is.na(percent_expressing), "n=0", "")),
            size = 2.6, color = "grey40", vjust = 1.8) +
  facet_grid(. ~ group_definition, scales = "free", space = "free_x") +
  scale_fill_viridis_c(name = "Mean expr\n(log-normalized)", option = "C",
                       direction = -1, na.value = "grey90") +
  scale_color_viridis_c(guide = "none", option = "C", direction = -1, na.value = "grey90") +
  scale_size_continuous(name = "% expressing", range = c(1.5, 8),
                        limits = c(0, 100)) +
  scale_x_discrete(labels = function(x) paste0(x, "\n(", sapply(x, function(s)
    as.character(sample_map[[s]])), ")")) +
  labs(x = NULL, y = NULL,
       title = NULL) +
  theme_classic(base_family = "Liberation Sans") +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

save_fig(p, "R2_Fig_canonical_markers_per_sample", w = 9, h = 3.8)

writeLines(capture.output(sessionInfo()),
           file.path(OUT_TAB, "R2_canonical_markers_per_sample_sessionInfo.txt"))
cat("\nRun finished:", format(Sys.time()), "\n")
cat("=== R2 Task 2 (round 2) Complete ===\n")
