# ==============================================================================
# R2_task1b_sensitivity_condition_composition.R
# Reviewer-2 integrity audit, Task 1b: minimal sensitivity analysis for the
# condition-composition imbalance between doublet comparators.
#
# Motivation (from R2_task1_cell_annotation_provenance.R):
#   A (Fibrocartilage-dominant) is 99.78% ONFH cells;
#   B (Other cartilage-lineage) is 99.62% Control cells.
#   The cell-level A-vs-B doublet comparison is therefore almost fully
#   confounded with disease condition. Donor is the biological unit; the
#   donor-level paired test (p=0.844) already showed no significant difference.
#
# Minimal sensitivity analyses (descriptive, no new claims):
#   S1: A vs C (Other fibrocartilage-candidate; 93.1% ONFH) - condition-balanced
#       program comparison: cell-level Fisher + donor-level paired Wilcoxon.
#   S2: A vs B restricted to ONFH donors (B has only 114 ONFH cells).
#   S3: doublet rate by condition within each lineage_call stratum
#       (quantifies the condition effect itself).
#
# Output: tables/R2_dblsens_condition_sensitivity.csv
# Seed 42. Previous results untouched.
# ==============================================================================

set.seed(42)
suppressMessages(library(Seurat))

db <- readRDS("/mnt/shared-workspace/shared/phase1/GSE311084_scored_doublet.rds")
md <- db@meta.data
md$dbl <- tolower(md$scDblFinder_class) == "doublet"
TAB <- "/mnt/results/phase1/tables"

A <- md$dominant_program == "Fibrocartilage"
Bl <- md$lineage_call == "cartilage_lineage" & !A
C <- md$lineage_call == "cartilage_fibrocartilage_candidate" & !A
onfh <- md$condition == "ONFH"

fisher_or <- function(sel1, sel2) {
  tab <- matrix(c(sum(md$dbl[sel1]), sum(!md$dbl[sel1]),
                  sum(md$dbl[sel2]), sum(!md$dbl[sel2])), nrow = 2,
                dimnames = list(c("doublet","singlet"), c("G1","G2")))
  ft <- fisher.test(tab)
  c(rate1 = 100 * mean(md$dbl[sel1]), n1 = sum(sel1),
    rate2 = 100 * mean(md$dbl[sel2]), n2 = sum(sel2),
    p = ft$p.value, OR = unname(ft$estimate))
}

# donor-level paired Wilcoxon for A vs C
pd <- do.call(rbind, lapply(sort(unique(md$donor_id)), function(d) {
  da <- md$dbl[A & md$donor_id == d]; dc <- md$dbl[C & md$donor_id == d]
  data.frame(donor = d, condition = unique(md$condition[md$donor_id == d]),
             n_A = length(da), pct_A = 100 * mean(da),
             n_C = length(dc), pct_C = 100 * mean(dc))
}))
pd$valid <- pd$n_A > 0 & pd$n_C > 0
wt <- suppressWarnings(wilcox.test(pd$pct_A[pd$valid], pd$pct_C[pd$valid], paired = TRUE))

rows <- list()
r <- fisher_or(A, C)
rows[[1]] <- data.frame(analysis = "S1_A_vs_C_cell_level_fisher", value = paste0(
  "A: ", sprintf("%.2f", r[["rate1"]]), "% (", r[["n1"]], " cells) vs C: ",
  sprintf("%.2f", r[["rate2"]]), "% (", r[["n2"]], " cells); Fisher p=",
  format(r[["p"]], digits = 3), ", OR=", sprintf("%.3f", r[["OR"]]),
  " [descriptive; cells not independent]"))
rows[[2]] <- data.frame(analysis = "S1_A_vs_C_donor_paired_wilcoxon", value = paste0(
  "n pairs = ", sum(pd$valid), "/", nrow(pd),
  "; p = ", format(wt$p.value, digits = 3),
  "; per-donor A%: ", paste(sprintf("%.2f", pd$pct_A), collapse = "/"),
  "; C%: ", paste(sprintf("%.2f", pd$pct_C), collapse = "/")))
r2 <- fisher_or(A & onfh, Bl & onfh)
rows[[3]] <- data.frame(analysis = "S2_A_vs_B_within_ONFH_cell_level_fisher", value = paste0(
  "A(ONFH): ", sprintf("%.2f", r2[["rate1"]]), "% (", r2[["n1"]],
  " cells) vs B(ONFH): ", sprintf("%.2f", r2[["rate2"]]), "% (", r2[["n2"]],
  " cells); Fisher p=", format(r2[["p"]], digits = 3), ", OR=",
  sprintf("%.3f", r2[["OR"]]), " [B(ONFH) n very small; unstable]"))
# S3: condition effect within strata
for (st in unique(c("cartilage_lineage", "cartilage_fibrocartilage_candidate"))) {
  sel <- md$lineage_call == st
  r3 <- fisher_or(sel & onfh, sel & !onfh)
  rows[[length(rows) + 1]] <- data.frame(analysis = paste0("S3_condition_within_", st), value = paste0(
    "ONFH: ", sprintf("%.2f", r3[["rate1"]]), "% (", r3[["n1"]],
    " cells) vs Control: ", sprintf("%.2f", r3[["rate2"]]), "% (", r3[["n2"]],
    " cells); Fisher p=", format(r3[["p"]], digits = 3), ", OR=",
    sprintf("%.3f", r3[["OR"]]), " [descriptive]"))
}
out <- do.call(rbind, rows)
print(out)
write.csv(out, file.path(TAB, "R2_dblsens_condition_sensitivity.csv"), row.names = FALSE)
cat("\nDONE task1b sensitivity script\n")
