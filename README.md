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

**Provenance note (honest disclosure):** Phase 1 was executed in an earlier working session. Session-local copies of scripts P1-00 through P1-06 and P1-08 through P1-10 were lost in a workspace storage reset; only `07_freeze_signatures.R` (preserved under `/mnt/results/`) survives verbatim. According to the original analysis report, Phase 1 results were persisted in the working environment (including Seurat objects, frozen signatures, tables and figures); this public repository snapshot does not contain all those inputs or full Seurat objects, and every parameter is documented in `parameters.md` and `FINAL_PHASE1_REPORT.md`. Many revision-phase scripts were executed in the original analysis environment; check the included run logs. Not all reconstructed Phase-1 modules are executable end-to-end.

**Reconstructed Phase-1 scripts (`scripts/reconstructed_phase1/`, added 2026-10-15):** partially reconstructed scripts or documented workflow skeletons for P1-00/01, P1-02/03, P1-04, P1-05, P1-06, P1-08, P1-09 and P1-10 were prepared from the persisted objects, documented parameters and saved result tables. Every file is headed `RECONSTRUCTED` and is **not** claimed to be the original script. Parameters that could not be recovered are marked `UNKNOWN` inside the files. Executed consistency checks (`RECONSTRUCTED_validation_check.R`, log in `logs/`):

| Check | Module | Result |
|---|---|---|
| V1 | P1-05 scoring | AddModuleScore recomputed from frozen curated signatures reproduces the stored `MS_1..MS_5` exactly (per-cell Spearman = 1.000 for all 5 programs; max donor-mean diff = 0.0000 vs `donor_module_scores.csv`) |
| V2 | P1-06 pseudobulk DE | edgeR recompute vs saved table: 18,450 genes, Spearman rho(logFC) = 0.9999; FDR<0.05 overlap Jaccard = 0.826 (2,340 vs 1,970 significant genes; consistent with edgeR version drift — logFC and effect direction fully reproduced) |
| V3 | P1-02/03 clustering | res 0.4 = 13 clusters as documented; single RNA assay confirms no batch integration |
| V4 | P1-00/01 QC | 57,476 QC'd − 166 phase-0.5 Macrophage-cluster cells = 57,310 final; per-donor exact (see `R2_cell_identity_audit.md`) |
| V5 | P1-05 AUCell | `AU_*` columns unpopulated (0/57,310), consistent with "AUCell omitted" |

Previously executed validations (revision phase): GSE74089 DE re-run rho(logFC) = 1.0 (`R2_task5_gse74089_rerun.R`); ssGSEA disease-effect reproduction rho = 0.982 (GSE74089) / 0.787 (GSE345158) (`R2_ext_specificity.R`). Steps without an executed re-run (raw QC filtering, integration-free clustering, lineage re-derivation, robustness re-execution) are marked **UNVERIFIED** in the reconstructed scripts; the persisted objects and tables remain the authoritative results.

**Annotation provenance (audit 2026-10-15):** the 57,310 GSE311084 cells are cartilage-lineage cells retained after the phase-0.5 marker-based lineage re-annotation of an integrated object (phase-0 per-sample preliminary annotations were heterogeneous: Fibroblast 38,585 / Chondrocyte 18,765 / Macrophage 126 of 57,476 QC'd cells; the 166 cells in the phase-0.5 Macrophage cluster were removed). See `R2_cell_identity_audit.md`.

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
| `R2_ext_specificity.R` | R2 Task 4 (ext.) | External-cohort specificity: ssGSEA of frozen Fibrocartilage signature vs 1,000 expression-matched random sets + EMT/ECM/fibroblast comparators in GSE74089 & GSE345158 | `R2_extspec_*.csv`, `R2_Fig_extspec.svg` |
| `R2_doublet_sensitivity.R` | R2 Task 3 (sens.) | Doublet sensitivity: doublet fraction in Fibrocartilage-dominant vs other cells; donor-level scores before/after doublet removal (MS_2 primary, UCell secondary); exact + Monte-Carlo permutation | `R2_dblsens_*.csv`, `R2_Fig_dblsens.svg` |
| `phase1_07_freeze_signatures.R` | (Phase 1) | Signature freezing prior to validation (preserved original) | `frozen_signatures.csv` |

## 3. Data Availability

| Dataset | Role | Source |
|---------|------|--------|
| GSE311084 | Discovery (scRNA-seq, 6 donors, 57,310 cartilage-lineage cells after phase-0.5 marker-based re-annotation) | GEO |
| GSE290411 | Cross-context (scRNA-seq, 4 SIONFH donors, 2,113 cells) | GEO |
| GSE74089 | External validation (Agilent microarray, hip cartilage, 4 vs 4) | GEO |
| GSE345158 | External validation (bulk RNA-seq, femoral head, 5 vs 4) | GEO |

All four datasets are public; no new human material was generated. Frozen signatures (215 genes, 8 signatures) are in `tables/frozen_signatures.csv` and must not be altered in any re-analysis.

## 4. Statistical Principles

- Biological unit = donor. GSE311084 is n=3 vs 3; cells are never treated as independent replicates.
- Primary inferential contrasts use donors as biological replicates (Welch t-test, Cohen's d, donor-level pseudobulk). Certain cell-level proportion tests are explicitly descriptive and do not constitute biological replication.
- External-validation disease effects combined by random-effects meta-analysis (DerSimonian-Laird; heterogeneity reported as I²).
- Signature disease effects tested against 1,000 expression-matched random gene sets (two-sided empirical p).
- No machine-learning classifiers, no drug-prediction, no network-pharmacology analyses are included.


## Release-specific reproducibility limitations (v2.0-revision-r2)

- This repository combines preserved original code, executed revision scripts, partial Phase-1 reconstructions, and supporting output tables. It is **not** a verified end-to-end recreation of the original Phase-1 processing.
- Some `RECONSTRUCTED_` scripts contain placeholder functions (`NULL`) and documentary pseudocode; they are not complete pipelines. QC thresholds and several original processing decisions are marked `UNKNOWN` where unavailable.
- Reconstructed edgeR log-fold changes were highly concordant but the FDR-significant set was not identical (Jaccard 0.826).
- Executed scripts depend on study inputs / RDS objects not shipped in this archive and use environment-specific absolute paths such as `/mnt/shared-workspace/shared` or `/workspace`. Adapt these paths after retrieving the inputs from GEO or the authors.
- The annotated local git commit `c7cf65b` was created by Biomni in its own workspace. Copying these files into an existing remote repository creates a different commit SHA; the new public version tag should be created at the actual pushed repository commit, not assumed to be the local Biomni SHA.
- Re-analysis should not silently change the pre-frozen gene signature.
