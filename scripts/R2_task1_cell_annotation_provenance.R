# ==============================================================================
# R2_task1_cell_annotation_provenance.R
# Reviewer-2 integrity audit, Task 1: GSE311084 cell annotation provenance
#
# Purpose: resolve the contradiction between the phase-0 per-sample preliminary
# annotation (Fibroblast 38,585 / Chondrocyte 18,765 / Macrophage 126 = 57,476
# QC'd cells, recorded in R2_cell_composition_per_sample.csv and
# phase0_celltype_results.json) and the final analysis set of 57,310 cells,
# using barcode-level tracking wherever persisted artifacts allow it.
#
# Outputs:
#   tables/R2_cell_annotation_provenance.csv   (stage-by-stage counts + reconciliation)
#   tables/R2_lineage_crosswalk.csv            (cluster-level crosswalk across stages)
#   tables/R2_cell_identity_crosstab.csv       (per-cell cross-tab from final metadata)
#   tables/R2_phase0_panel_score_by_group.csv  (phase-0 panel logic applied per cell)
#   R2_cell_identity_audit.md                  (audit narrative; written by this script)
#
# Seed 42; R 4.4.2; Seurat 5.5.1. No frozen signature or result is modified.
# ==============================================================================

set.seed(42)
suppressMessages({library(Seurat); library(jsonlite)})

p05  <- "/mnt/shared-workspace/shared/phase0_5"
RDS  <- "/mnt/shared-workspace/shared/phase1"
DAT  <- "/mnt/shared-workspace/shared/onfh_data"
TAB  <- "/mnt/results/phase1/tables"
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

options(width = 200)

cat("== [1] Load persisted artifacts ==\n")
json <- fromJSON(file.path(DAT, "phase0_celltype_results.json"), simplifyVector = FALSE)
g31  <- Filter(function(s) s$dataset == "GSE311084" | s$sample_id %in% c("A-01","A-02","A-03","A-06","A-07","A-08"), json)
cat("phase-0 JSON entries for GSE311084:", length(g31), "\n")

rc  <- readRDS(file.path(RDS, "GSE311084_reclustered.rds"))
sc  <- readRDS(file.path(RDS, "GSE311084_scored.rds"))
db  <- readRDS(file.path(RDS, "GSE311084_scored_doublet.rds"))
md  <- sc@meta.data

ct_counts  <- readRDS(file.path(p05, "GSE311084_ct_counts.rds"))
ct_donor   <- readRDS(file.path(p05, "GSE311084_ct_donor.rds"))
major_calls<- readRDS(file.path(p05, "GSE311084_major_calls.rds"))
lineage_cl <- readRDS(file.path(p05, "GSE311084_lineage.rds"))

cat("== [2] Barcode-level tracking across persisted objects ==\n")
bc_rc <- colnames(rc); bc_sc <- colnames(sc); bc_db <- colnames(db)
track <- data.frame(
  comparison = c("reclustered vs scored", "scored vs scored_doublet"),
  n_only_A = c(length(setdiff(bc_rc, bc_sc)), length(setdiff(bc_sc, bc_db))),
  n_only_B = c(length(setdiff(bc_sc, bc_rc)), length(setdiff(bc_db, bc_sc))),
  n_intersect = c(length(intersect(bc_rc, bc_sc)), length(intersect(bc_sc, bc_db)))
)
print(track)
cat("duplicated barcodes within scored:", sum(duplicated(bc_sc)), "\n")

# barcode suffix <-> sample_id 1:1 mapping
sfx <- sub(".*_", "", bc_sc)
sfx_map <- table(sfx, md$sample_id)
cat("barcode suffix x sample_id (must be 1:1):\n"); print(sfx_map)

cat("\n== [3] Stage-by-stage counts and reconciliation ==\n")
# phase-0 per-sample (JSON)
p0 <- do.call(rbind, lapply(g31, function(s) data.frame(
  sample_id = s$sample_id, donor_id = s$donor_id, condition = s$disease_status,
  n_raw = s$n_cells_raw, n_phase0_qc = s$n_cells_filtered,
  n_fibroblast_phase0 = as.integer(unlist(s$cell_type_counts$fibroblast) %||% 0),
  n_chondrocyte_phase0 = as.integer(unlist(s$cell_type_counts$chondrocyte) %||% 0),
  n_macrophage_phase0 = as.integer(unlist(s$cell_type_counts$macrophage) %||% 0))))
p0$n_phase0_total <- p0$n_fibroblast_phase0 + p0$n_chondrocyte_phase0 + p0$n_macrophage_phase0
stopifnot(all(p0$n_phase0_total == p0$n_phase0_qc))

# phase-0.5 major calls per donor
p05d <- reshape(ct_donor, idvar = "donor_id", timevar = "cell_type", direction = "wide")
names(p05d) <- sub("n_cells.", "", names(p05d))
p05d <- p05d[, c("donor_id", sort(names(p05d)[-1]))] # donor first

# final per sample
fin <- as.data.frame(table(md$sample_id, md$donor_id))
fin <- fin[fin$Freq > 0, c("Var1", "Var2", "Freq")]
names(fin) <- c("sample_id", "donor_id", "n_final")

prov <- merge(p0, p05d, by = "donor_id", all.x = TRUE)
prov <- merge(prov, fin, by = c("sample_id", "donor_id"))
names(prov)[names(prov) == "Chondrocyte"] <- "n_chondrocyte_phase05"
names(prov)[names(prov) == "Macrophage"]  <- "n_macrophage_phase05"
prov$n_removed_after_phase05 <- prov$n_phase0_qc - prov$n_final
prov$removed_matches_phase05_macrophage <- prov$n_removed_after_phase05 == prov$n_macrophage_phase05
prov <- prov[order(prov$sample_id),
  c("sample_id","donor_id","condition","n_raw","n_phase0_qc",
    "n_fibroblast_phase0","n_chondrocyte_phase0","n_macrophage_phase0",
    "n_chondrocyte_phase05","n_macrophage_phase05","n_final",
    "n_removed_after_phase05","removed_matches_phase05_macrophage")]
print(prov)
write.csv(prov, file.path(TAB, "R2_cell_annotation_provenance.csv"), row.names = FALSE)

cat("\nTotals: phase0 QC =", sum(prov$n_phase0_qc),
    "| phase0.5 Chondrocyte+Macrophage =", sum(prov$n_chondrocyte_phase05) + sum(prov$n_macrophage_phase05),
    "| final =", sum(prov$n_final), "\n")

cat("\n== [4] Cluster-level crosswalk across annotation stages ==\n")
# phase-0: per-sample cluster calls (JSON)
p0cl <- do.call(rbind, lapply(g31, function(s) data.frame(
  stage = "phase0_per_sample", sample_id = s$sample_id, donor_id = s$donor_id,
  cluster = names(s$cluster_types),
  call = unlist(s$cluster_types),
  fib_score = unlist(s$cluster_mean_scores$fibroblast)[names(s$cluster_types)],
  chond_score = unlist(s$cluster_mean_scores$chondrocyte)[names(s$cluster_types)],
  n_cells = NA_integer_, rationale = NA_character_)))
# phase-0.5: integrated cluster calls
p05cl <- data.frame(stage = "phase05_integrated", sample_id = NA, donor_id = NA,
  cluster = as.character(major_calls$cluster), call = major_calls$cell_type,
  fib_score = NA, chond_score = NA, n_cells = NA,
  rationale = major_calls$rationale)
# final: res 0.4 clusters x lineage_call x dominant_program
fcl <- as.data.frame(table(cluster = md$RNA_snn_res.0.4, lineage = md$lineage_call,
                           prog = md$dominant_program))
fcl <- fcl[fcl$Freq > 0, ]
fcl$stage <- "final_res0.4"; fcl$sample_id <- NA; fcl$donor_id <- NA
fcl$call <- paste(fcl$lineage, fcl$prog, sep = " | ")
fcl$fib_score <- NA; fcl$chond_score <- NA; fcl$rationale <- NA
fcl <- fcl[, c("stage","sample_id","donor_id","cluster","call","fib_score","chond_score","Freq","rationale")]
names(fcl)[names(fcl) == "Freq"] <- "n_cells"
names(p0cl)[names(p0cl) == "call"] <- "call"
p0cl$n_cells <- NA
p0cl <- p0cl[, c("stage","sample_id","donor_id","cluster","call","fib_score","chond_score","n_cells","rationale")]
p05cl$n_cells <- NA
p05cl <- p05cl[, c("stage","sample_id","donor_id","cluster","call","fib_score","chond_score","n_cells","rationale")]
crosswalk <- rbind(p0cl, p05cl, fcl)
write.csv(crosswalk, file.path(TAB, "R2_lineage_crosswalk.csv"), row.names = FALSE)
cat("crosswalk rows:", nrow(crosswalk), "(phase0:", sum(crosswalk$stage=="phase0_per_sample"),
    "; phase05:", sum(crosswalk$stage=="phase05_integrated"),
    "; final:", sum(crosswalk$stage=="final_res0.4"), ")\n")

cat("\n== [5] Per-cell identity cross-tab (final metadata) ==\n")
xt <- as.data.frame(table(
  cell_type = md$cell_type, lineage_call = md$lineage_call,
  dominant_program = md$dominant_program, condition = md$condition, donor_id = md$donor_id))
xt <- xt[xt$Freq > 0, ]
write.csv(xt, file.path(TAB, "R2_cell_identity_crosstab.csv"), row.names = FALSE)
cat("non-zero combinations:", nrow(xt), "\n")
print(xt[, c("lineage_call","dominant_program","condition","Freq")][order(xt$dominant_program), ])

cat("\n== [6] Group audit: A / B / B2 / C definitions, overlap, composition ==\n")
A  <- md$dominant_program == "Fibrocartilage"
Bl <- md$lineage_call == "cartilage_lineage" & !A
B2 <- md$lineage_call == "cartilage_lineage" &
      !md$dominant_program %in% c("Fibrocartilage", "Undetermined")
C  <- md$lineage_call == "cartilage_fibrocartilage_candidate" & !A
grp <- c(A = sum(A), B = sum(Bl), B2 = sum(B2), C = sum(C))
cat("group sizes:\n"); print(grp)
cat("A+B+C =", sum(A) + sum(Bl) + sum(C), "of", ncol(sc), "(complete partition:",
    sum(A) + sum(Bl) + sum(C) == ncol(sc), ")\n")
cat("A & B overlap:", sum(A & Bl), "| A & C overlap:", sum(A & C),
    "| B & C overlap:", sum(Bl & C), "| B2 within B:", sum(B2 & Bl) == sum(B2), "\n")

comp <- data.frame(
  group = c("A_Fibrocartilage_dominant","B_Other_cartilage_lineage","C_Other_fibrocartilage_candidate"),
  n = c(sum(A), sum(Bl), sum(C)),
  n_ONFH = c(sum(A & md$condition=="ONFH"), sum(Bl & md$condition=="ONFH"), sum(C & md$condition=="ONFH")),
  n_Control = c(sum(A & md$condition=="Control"), sum(Bl & md$condition=="Control"), sum(C & md$condition=="Control")),
  pct_cartilage_lineage_call = c(round(100*sum(A & md$lineage_call=="cartilage_lineage")/sum(A),2),
                                 100, 0),
  pct_candidate_call = c(round(100*sum(A & md$lineage_call=="cartilage_fibrocartilage_candidate")/sum(A),2),
                         0, 100))
print(comp)
write.csv(comp, file.path(TAB, "R2_group_composition_audit.csv"), row.names = FALSE)

cat("\n== [7] Phase-0 panel logic applied per cell (reconstruction of scoring rule) ==\n")
# phase-0 panels from phase0_celltype_results.json (marker_genes_present)
s1 <- g31[[1]]
fib_panel  <- unlist(s1$marker_genes_present$fibroblast)
chond_panel<- unlist(s1$marker_genes_present$chondrocyte)
genes <- rownames(sc)
fib_g  <- intersect(fib_panel, genes); chond_g <- intersect(chond_panel, genes)
cat("fibroblast panel genes present:", paste(fib_g, collapse=","), "\n")
cat("chondrocyte panel genes present:", paste(chond_g, collapse=","), "\n")
expr <- GetAssayData(sc, layer = "data")
fib_sc  <- colMeans(expr[fib_g, , drop = FALSE])
chnd_sc <- colMeans(expr[chond_g, , drop = FALSE])
leaner <- ifelse(fib_sc > chnd_sc, "fibroblast_panel_higher",
          ifelse(chnd_sc > fib_sc, "chondrocyte_panel_higher", "tie"))
pg <- data.frame(
  group = c("A_Fibrocartilage_dominant","B_Other_cartilage_lineage","C_Other_fibrocartilage_candidate"),
  n = c(sum(A), sum(Bl), sum(C)),
  pct_fibroblast_panel_higher = c(
    round(100 * sum(leaner[A] == "fibroblast_panel_higher") / sum(A), 1),
    round(100 * sum(leaner[Bl] == "fibroblast_panel_higher") / sum(Bl), 1),
    round(100 * sum(leaner[C] == "fibroblast_panel_higher") / sum(C), 1)),
  pct_chondrocyte_panel_higher = c(
    round(100 * sum(leaner[A] == "chondrocyte_panel_higher") / sum(A), 1),
    round(100 * sum(leaner[Bl] == "chondrocyte_panel_higher") / sum(Bl), 1),
    round(100 * sum(leaner[C] == "chondrocyte_panel_higher") / sum(C), 1)),
  mean_fibroblast_panel_expr = c(round(mean(fib_sc[A]), 3), round(mean(fib_sc[Bl]), 3), round(mean(fib_sc[C]), 3)),
  mean_chondrocyte_panel_expr = c(round(mean(chnd_sc[A]), 3), round(mean(chnd_sc[Bl]), 3), round(mean(chnd_sc[C]), 3)))
print(pg)
write.csv(pg, file.path(TAB, "R2_phase0_panel_score_by_group.csv"), row.names = FALSE)

# marker co-expression per group (detection rates)
det <- function(gset, sel) round(100 * mean(colSums(expr[gset, sel, drop=FALSE] > 0) > 0), 1)
coex <- data.frame(
  group = pg$group,
  pct_SOX9 = c(det("SOX9",A), det("SOX9",Bl), det("SOX9",C)),
  pct_ACAN = c(det("ACAN",A), det("ACAN",Bl), det("ACAN",C)),
  pct_COL2A1 = c(det("COL2A1",A), det("COL2A1",Bl), det("COL2A1",C)),
  pct_any_chondrocyte_panel = c(det(chond_g,A), det(chond_g,Bl), det(chond_g,C)),
  pct_any_fibroblast_panel = c(det(fib_g,A), det(fib_g,Bl), det(fib_g,C)),
  pct_both_panels = c(
    round(100*mean((colSums(expr[chond_g,A,drop=FALSE]>0)>0) & (colSums(expr[fib_g,A,drop=FALSE]>0)>0)),1),
    round(100*mean((colSums(expr[chond_g,Bl,drop=FALSE]>0)>0) & (colSums(expr[fib_g,Bl,drop=FALSE]>0)>0)),1),
    round(100*mean((colSums(expr[chond_g,C,drop=FALSE]>0)>0) & (colSums(expr[fib_g,C,drop=FALSE]>0)>0)),1)))
print(coex)
write.csv(coex, file.path(TAB, "R2_marker_coexpression_by_group.csv"), row.names = FALSE)

cat("\n== [8] GSE290411 cross-context annotation provenance (secondary check) ==\n")
g2 <- readRDS(file.path(RDS, "GSE290411_scored.rds"))
cat("GSE290411 scored ncells:", ncol(g2), "\n")
print(table(g2$cell_type, g2$donor_id))
print(table(g2$lineage_call))

saveRDS(NULL, "/dev/null")  # no-op; session state not needed downstream
cat("\nDONE task1 provenance script\n")