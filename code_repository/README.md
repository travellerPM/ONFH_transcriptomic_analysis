# Code Repository: Cross-Cohort Single-Cell Transcriptomics of ONFH Chondrocyte-State Remodeling

Manuscript: Conserved and context-specific chondrocyte-state remodeling in osteonecrosis of the femoral head (ONFH) — cross-cohort single-cell and bulk transcriptomic analysis.

Revision phase code package (Reviewer 2, Task 6: code and parameter availability).

## 1. Environment

- R 4.4.2 (see `sessionInfo.txt` for full package versions)
- Key packages: Seurat 5.5.1, UCell, scDblFinder, edgeR, DESeq2, limma, fgsea, decoupleR, dorothea, metafor, openxlsx
- All scripts assume `R_LIBS_USER=/workspace/.Rlib` and data under `/mnt/shared-workspace/shared/`

## 2. Pipeline Overview and Run Order

### Phase 1 — Discovery, cross-context, validation, robustness (original submission)

| Step | Input | Output (persisted) | Key parameters |
|------|-------|--------------------|----------------|
| P1-00/01 Raw data + QC | GSE311084, GSE290411 raw counts (GEO) | `phase0_celltype_results.json`, Seurat objects | Per-sample QC; filtered cell counts recorded per sample |
| P1-02/03 Normalization, clustering | Filtered objects | `*_reclustered.rds` | SCT/LogNorm; GSE311084 res=0.4 (13 clusters); GSE290411 res=0.3 (6 clusters); see `resolution_selection_summary.md` |
| P1-04 Annotation | Clustered objects | `lineage_call`, `dominant_program` metadata | Marker-based; cartilage_lineage vs cartilage_fibrocartilage_candidate |
| P1-05 Program scoring | Annotated objects | `*_scored.rds`, `donor_module_scores.csv` | AddModuleScore + UCell (5 programs); AUCell omitted (memory, 57k cells) |
| P1-06 Pseudobulk DE | Scored object | `pseudobulk_DE_edgeR.csv` | edgeR primary (donor-level pseudobulk, n=3 vs 3); DESeq2 + limma-voom concordance |
| P1-07 Signature freezing | Scored object + edgeR | `frozen_signatures.csv`, `frozen_chondrocyte_signatures.xlsx` | **Script preserved: `scripts/phase1_07_freeze_signatures.R`** — frozen BEFORE examining validation data |
| P1-08 External validation | GSE74089, GSE345158 | `validation_grading_summary.csv` | GSE74089: limma + ssGSEA + mean z-score; GSE345158: DESeq2 + ssGSEA + mean z-score |
| P1-09 Robustness | Scored objects | robustness tables | LODO (6 iterations), downsampling 25–90%, AddModuleScore vs UCell, leave-one-gene-out, 7 resolutions |
| P1-10 Figures | All tables | `figures/Fig1–Fig7` | ggplot2, SVG + PNG 300 dpi |

**Provenance note (honest disclosure):** Phase 1 was executed in an earlier working session. Session-local copies of scripts P1-00 through P1-06 and P1-08 through P1-10 were lost in a workspace storage reset; only `07_freeze_signatures.R` (preserved under `/mnt/results/`) survives verbatim. All Phase 1 *outputs* are fully persisted (Seurat objects, frozen signatures, 54 result tables, 14 figures), and every parameter is documented in `parameters.md` and `FINAL_PHASE1_REPORT.md`. The revision-phase scripts below are complete and were all executed in the current session.

### Phase 2 — Revision analyses (current package, all scripts executable)

Run order matches reviewer-response numbering:

| Script | Reviewer point | Purpose | Outputs |
|--------|----------------|---------|---------|
| `R2_task3_doublet.R` | R2 Task 3 | scDblFinder doublet detection, both sc datasets, per-donor | `R2_doublet_rates_*.csv`, `R2_doublet_QC_*.csv`, per-cell calls, `*_scored_doublet.rds` |
| `R2_task4_negative_control.R` | R2 Task 4 | 1,000 expression-matched random gene sets per program; empirical p-values; comparison signatures (EMT, generic ECM, generic fibroblast) | `R2_negative_control_results.csv`, null distributions, `R2_Fig_negative_control.svg` |
| `R1_meta_analysis.R` | R1 #2 | Random-effects (DerSimonian-Laird) meta-analysis of disease effect across GSE74089 + GSE345158; forest plot | `R1_meta_analysis_results.csv`, `R1_Fig_meta_forest.svg` |
| `R2_task1_donor_figure.R` | R2 Task 1 | Donor-level module-score figure + Welch t-test / Cohen's d at donor level (n=3 vs 3) | `R2_donor_scores_*_stats.csv`, `R2_Fig_donor_scores.svg` |
| `R2_task2_cell_composition.R` | R2 Task 2 | Per-sample cell composition; canonical marker confirmation of chondrocyte identity | `R2_cell_composition_per_sample.csv`, `R2_canonical_markers_summary.csv` |
| `R2_task5_gse74089_rerun.R` | R2 Task 5 | Independent re-run of GSE74089 probe annotation + DE; verification against original | `R2_gse74089_rerun_verification.csv` (Spearman rho = 1.0) |
| `R1_reproducibility_framework.R` | R1 #1 | Consolidated 12-item reproducibility evidence table | `R1_reproducibility_framework_summary.csv` |
| `R1_tissue_context.R` | R1 #3 | Tissue-context interpretation of validation cohorts | `R1_tissue_context_*.csv`, `R1_tissue_context_interpretation.xlsx` |
| `phase1_07_freeze_signatures.R` | (Phase 1) | Signature freezing prior to validation (preserved original) | `frozen_signatures.csv` |

## 3. Data Availability

| Dataset | Role | Source |
|---------|------|--------|
| GSE311084 | Discovery (scRNA-seq, 6 donors, 57,310 chondrocytes) | GEO |
| GSE290411 | Cross-context (scRNA-seq, 4 SIONFH donors, 2,113 cells) | GEO |
| GSE74089 | External validation (Agilent microarray, hip cartilage, 4 vs 4) | GEO |
| GSE345158 | External validation (bulk RNA-seq, femoral head, 5 vs 4) | GEO |

All four datasets are public; no new human material was generated. Frozen signatures (215 genes, 8 signatures) are in `tables/frozen_signatures.csv` and must not be altered in any re-analysis.

## 4. Statistical Principles

- Biological unit = donor. GSE311084 is n=3 vs 3; cells are never treated as independent replicates.
- All group comparisons at donor level (Welch t-test, Cohen's d) or donor-level pseudobulk (edgeR/DESeq2/limma-voom).
- External-validation disease effects combined by random-effects meta-analysis (DerSimonian-Laird; heterogeneity reported as I²).
- Signature disease effects tested against 1,000 expression-matched random gene sets (two-sided empirical p).
- No machine-learning classifiers, no drug-prediction, no network-pharmacology analyses are included.
