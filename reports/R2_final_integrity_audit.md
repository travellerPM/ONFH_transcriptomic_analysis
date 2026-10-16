# R2 Final Integrity Audit — Reviewer-2 Revision Analyses

**Date:** 2026-10-15
**Scope:** evidence-chain gaps, annotation contradictions, and code reproducibility in the existing ONFH revision analyses. No frozen signature, sample set, or primary hypothesis was changed. All statistics use the donor as the biological unit.

---

## Task 1 — GSE311084 cell annotation provenance (highest priority)

**Original problem.** `R2_cell_composition_per_sample.csv` records phase-0 per-sample preliminary annotations (Fibroblast 38,585 / Chondrocyte 18,765 / Macrophage 126 = 57,476 QC'd cells), while `R2_task3_doublet_cartilage_lineage_comparator.R` claimed "All cells in GSE311084 are phase-0 Chondrocyte calls" for the 57,310-cell final set. The two statements could not both be correct.

**Verification method.** Barcode-level set operations across all persisted objects; parsing of `phase0_celltype_results.json` (per-sample cluster-level calls and 13-lineage marker panels); inspection of phase-0.5 artifacts (`GSE311084_major_calls.rds`, `GSE311084_lineage.rds`, `GSE311084_ct_counts.rds`, `GSE311084_ct_donor.rds`); per-cell cross-tabulation of the final metadata; rule-based reconstruction of the phase-0 scoring logic per cell.

**Actual findings.**
1. The two descriptions refer to different stages. Phase-0 = per-sample preliminary annotation (Fibroblast-majority). Phase-0.5 = integrated re-annotation with a 9-lineage marker panel: clusters 0–9 → Chondrocyte (`cartilage_lineage` 29,915 / `cartilage_fibrocartilage_candidate` 27,395 per cell), cluster 10 → Macrophage (166 cells; ONFH-1: 7, ONFH-2: 159). The 166 Macrophage-cluster cells were removed before the final object: 57,476 − 166 = 57,310, reconciled exactly per donor (6/6).
2. **The original claim was wrong as stated and has been corrected** (dated correction notes in the script and in `R2_task2_task3_completion_report.md`): the correct statement is "all 57,310 retained cells are Chondrocyte calls from the phase-0.5 integrated re-annotation". Terminology was not changed to mask the contradiction; the stage structure is documented.
3. Barcode-level verification: `reclustered` = `scored` = `scored_doublet` (57,310 barcodes, 0 differences, no duplicates); barcode suffix ↔ sample 1:1; donor/sample/condition mapping matches the phase-0 JSON.
4. Not recoverable (stated honestly): per-cell phase-0 labels (per-sample phase-0 objects lost; JSON is cluster-level) and the 166 removed barcodes (phase-0.5 intermediate object not persisted). Their removal is verified by exact per-donor count reconciliation, not barcode setdiff.
5. Group audit: A (7,149) + B (29,899) + C (20,262) = 57,310 — complete partition, zero overlap; B2 ⊂ B. A is 99.78% ONFH cells and 99.78% `cartilage_fibrocartilage_candidate`; B is 99.62% Control cells and 100% `cartilage_lineage`.
6. Phase-0 panel logic per cell: 100% of A-group cells score higher on the phase-0 fibroblast panel (COL1A1/2, COL3A1, DCN, LUM, PDGFRA) than the chondrocyte panel, yet 97.0% co-express both panels (SOX9 53.6%, ACAN 75.5%, COL2A1 60.5% detection). The phase-0 "fibroblast" majority was a panel artifact (no fibrocartilage category; DCN/LUM/COL1A1 are fibrocartilage ECM genes).
7. Minimal sensitivity analysis for the A/B composition imbalance (`R2_task1b_sensitivity_condition_composition.R`): A vs C (condition-balanced) doublet comparison 8.09% vs 8.91% (Fisher p = 0.034, OR = 0.899; donor-level paired Wilcoxon p = 0.59); A vs B within ONFH donors reverses direction (OR = 9.9; B(ONFH) n = 114). The cell-level A-vs-B doublet "depletion" is condition-composition-driven, not program-driven.

**Fixed?** Yes — corrected claim, corrected script comment, sensitivity analysis added.
**Impact on conclusions.** None on the primary donor-level result (Control 0.556 ± 0.141 vs ONFH 1.830 ± 0.267, Welch p = 0.0051; unchanged by doublet removal). The manuscript should describe the cells as cartilage-lineage after marker-based re-annotation (not "phase-0 chondrocytes") and should rest the doublet-robustness argument on donor-level tests and the A-vs-C sensitivity.
**Remaining limitations.** Per-cell phase-0 labels and the 166 barcodes are unrecoverable; the phase-0.5 lineage-calling clustering (11 clusters) persists only as cluster-level calls; scDblFinder calls are predictions.
**Evidence.** `R2_cell_identity_audit.md`; `tables/R2_cell_annotation_provenance.csv`, `R2_lineage_crosswalk.csv`, `R2_cell_identity_crosstab.csv`, `R2_group_composition_audit.csv`, `R2_phase0_panel_score_by_group.csv`, `R2_marker_coexpression_by_group.csv`, `R2_dblsens_condition_sensitivity.csv`; `logs/R2_task1_cell_annotation_provenance.log`, `logs/R2_task1b_sensitivity_condition_composition.log`.

---

## Task 2 — Phase-1 code repository repair (highest priority)

**Original problem.** README documented that P1-00…P1-06 and P1-08…P1-10 scripts were lost in a workspace storage reset; only `07_freeze_signatures.R` survives. Reviewer 2 requires a complete reproducible pipeline.

**Verification method.** Exhaustive search of the workspace, results directories, shared storage, and logs (only `07_freeze_signatures.R` found). Reconstruction from persisted objects, documented parameters (`parameters.md`, `FINAL_PHASE1_REPORT.md`), and saved result tables.

**Actual findings / actions.**
1. Eight equivalent executable scripts rebuilt under `scripts/reconstructed_phase1/`, every file headed `RECONSTRUCTED` and none claimed to be the original; unrecoverable parameters marked `UNKNOWN` inside the files (QC thresholds, integration-free clustering details, phase-0.5 marker-gene lists, PC count).
2. Executed consistency checks (`RECONSTRUCTED_validation_check.R`, EXIT 0):
   - V1 (P1-05 scoring): AddModuleScore recomputed from the frozen curated signatures reproduces the stored `MS_1…MS_5` **exactly** (per-cell Spearman = 1.000 for all 5 programs; max donor-mean diff = 0.0000 vs `donor_module_scores.csv`).
   - V2 (P1-06 pseudobulk DE): edgeR recompute vs saved table, 18,450 genes: Spearman rho(logFC) = 0.9999; FDR<0.05 Jaccard = 0.826 (2,340 vs 1,970 significant genes; consistent with edgeR version drift — logFC and directions fully reproduced).
   - V3 (P1-02/03): res 0.4 = 13 clusters as documented; single RNA assay confirms no batch integration.
   - V4 (P1-00/01): QC count reconciliation exact (Task 1).
   - V5: `AU_*` columns unpopulated (0/57,310), consistent with "AUCell omitted".
   - Previously executed validations: GSE74089 DE re-run rho = 1.0 (`R2_task5`); ssGSEA disease-effect reproduction rho = 0.982 / 0.787 (`R2_ext_specificity`).
3. README updated (reconstruction table, validation status, annotation-provenance correction, "57,310 cartilage-lineage cells" wording).
4. Git: no repository or remote existed. A local repository was initialized for the analysis working copy (commit `c7cf65b`, annotated tag `v2.0-revision-r2`). **Nothing has been pushed or published online**; `RELEASE_NOTES_v2.0-revision-r2.md` provides the prepared push/tag commands for the authors.

**Fixed?** Yes (reconstruction + executed validation + README + local version control).
**Impact on conclusions.** None — persisted objects/tables remain authoritative; the frozen-signature script is the original.
**Remaining limitations.** Raw-QC filtering, clustering, lineage re-derivation, and robustness re-execution are UNVERIFIED (marked in the scripts); edgeR FDR counts show version drift.
**Evidence.** `scripts/reconstructed_phase1/*`, `tables/RECONSTRUCTED_validation_results.csv`, `logs/RECONSTRUCTED_validation_check.log`, updated `README.md`, `RELEASE_NOTES_v2.0-revision-r2.md`.

---

## Task 3 — External-cohort negative-control statistical review

**Original problem.** (3.1) The directional-concordance null used the up-fraction of random-set genes, while the observed statistic matches genes to the frozen direction pattern (22 up + DCN down among 23 present genes) — definitionally asymmetric. (3.2) Raw ssGSEA disease effects were compared across signatures of very different size (Fibrocartilage 23 genes vs EMT 200 / ECM 111 / Fibroblast 88).

**Verification method.** Deterministic regeneration of the 1,000 expression-matched random sets (seed 42) — reproduced the saved null distribution exactly (max abs diff 4.44e-16, both cohorts); recomputation of the observed/null statistics from saved DE tables and per-sample scores; no ssGSEA re-run; original results preserved side by side.

**Actual findings.**
1. The ssGSEA disease-effect empirical p is definitionally symmetric (same set-level statistic on observed and random sets); recomputed p identical to saved (GSE74089 0.0010; GSE345158 0.0020). No correction required.
2. Corrected symmetric concordance null (same 22-up/1-down direction pattern assigned to matched random sets): GSE74089 observed 0.957 vs null mean 0.577, corrected p = 0.0010; GSE345158 observed 0.739 vs null 0.387, corrected p = 0.0030 — **conclusions unchanged** (original p = 0.0020 / 0.0030).
3. Direction-permutation null (permuting the frozen direction labels within the signature): p = 1.0 (GSE74089) and 0.326 (GSE345158) — the specific frozen direction pattern adds no information beyond coordinated up-regulation; **DCN's frozen down-direction is not reproduced in either cohort** (all 23 present genes are up in GSE74089). The concordance signal should be described as "coordinated up-regulation of the signature genes in disease", not as reproduction of the full direction pattern.
4. Effect-size standardization: null-standardized z for Fibrocartilage = 5.37 (GSE74089) / 4.64 (GSE345158); Cohen's d — GSE74089: Fibrocartilage 6.58, EMT 9.47, ECM 14.39, Fibroblast 2.74; GSE345158: Fibrocartilage 0.90, EMT −0.22, ECM 0.13, Fibroblast 0.82. Raw ssGSEA effects are not comparable across set sizes.
5. Three significance types separated: (a) empirical significance vs expression-matched random sets (Fibrocartilage p = 0.001 / 0.002); (b) disease-vs-control t-test (GSE74089 p = 0.0025; **GSE345158 p = 0.253, not significant**); (c) relative specificity vs generic programs (GSE74089: ECM/EMT show larger disease effects than Fibrocartilage; GSE345158: Fibrocartilage has the largest positive d). Random-set empirical p-values must not be presented as independent case-control p-values.

**Fixed?** Yes — corrected statistics reported alongside originals (no silent overwrite); interpretable standardized metrics added.
**Impact on conclusions.** Empirical-specificity conclusions stand; the manuscript must not claim (i) reproduction of the frozen direction pattern (DCN) or (ii) case-control significance in GSE345158 from the random-set p-value.
**Remaining limitations.** ssGSEA scores are platform-normalized but sample sizes are small (4v4, 5v4); comparator nulls were not re-generated (standardization uses t/d/z instead).
**Evidence.** `tables/R2_extspec_stats_review.csv`, `R2_extspec_corrected_concordance.csv`; `logs/R2_task3_extspec_stats_review.log`; script `R2_task3_extspec_stats_review.R`.

---

## Task 4 — Random-effects synthesis audit

**Original problem.** Verify effect-size metric, SE computation, DL estimator, τ²/I², cross-cohort scale comparability, and the appropriateness of pooling raw mean differences with only two cohorts.

**Verification method.** Line-by-line audit of `R1_meta_analysis.R`; exact recomputation from saved inputs; Hedges-g SMD sensitivity; fixed-effect comparison; descriptive synthesis.

**Actual findings.**
1. Reproduction is exact: max |original − recomputed| DL combined effect = 0 across all 8 signatures.
2. Metric = raw mean difference of ssGSEA scores (within-sample rank-normalized, comparable [0,1] scale across platforms — defensible); SE = |diff/t| is algebraically exact; DL correctly called via `metafor::rma(method="DL")`.
3. **Material finding:** with k = 2 (df = 1), τ²/I² are unstable. For Fibrocartilage: DL combined = 0.197 [−0.042, 0.436], **p = 0.107 (not significant)**, I² = 93.1%, τ² = 0.028; the fixed-effect estimate (0.247, p = 9.3e-18) is invalid under this heterogeneity. Hedges g: 5.72 (GSE74089) vs 0.836 (GSE345158); SMD-DL combined g = 3.00 [−1.76, 7.75], p = 0.216. Inflammatory is direction-discordant (−0.315 vs +0.050).
4. Prudent descriptive synthesis (recommended for the manuscript): per-cohort effects with CIs, direction concordance (Fibrocartilage 2/2 positive), per-cohort p-values, and the DL estimate reported with explicit k = 2 instability caveats; no substantive interpretation of I²/τ²; no claim of a significant pooled effect.

**Fixed?** Yes — audit table + descriptive synthesis + Methods/Discussion text provided; original results preserved.
**Impact on conclusions.** The pooled random-effects test for Fibrocartilage is not significant; external-validation support rests on per-cohort concordance of direction and the GSE74089 effect, with GSE345158 non-significant at case-control level. This must be reflected in the Response and manuscript.
**Remaining limitations.** Two cohorts, small n; SMD sensitivity available only for signatures with saved per-sample scores.
**Evidence.** `tables/R1_meta_audit_results.csv`; `logs/R2_task4_meta_analysis_audit.log`; script `R2_task4_meta_analysis_audit.R`.

---

## Task 5 — Figure 2 scoring-system check

**Original problem.** Results 3.2 quotes 0.556 vs 1.830 (donor-level AddModuleScore MS_2), while Figure 2 panels use UCell; the revision donor figure plotted AddModuleScore data labeled "UCell score".

**Verification method.** SVG text/label extraction and value-range checks against donor-level tables for every score figure; visual inspection; data-source audit of `R2_task1_donor_figure.R`.

**Actual findings.**
1. `Fig2_module_scores.svg` and `Fig2b_umap_feature_scores.svg`: cell-level UCell (title and 0–1 non-negative scale) — internally correctly labeled, but inconsistent with the AddModuleScore statistics quoted in the text.
2. `R2_Fig_donor_scores*.svg`: data are AddModuleScore donor means (match the text) but labeled "Donor mean UCell score" — **label bug**, corrected.
3. Re-outputs (originals preserved): `R2_Fig_donor_scores_corrected.*`, `R2_Fig_donor_scores_bar_corrected.*` (label fix, data unchanged); `Fig2_module_scores_donor_AddModuleScore.*` (donor-level MS_2, matches the text: 0.556 vs 1.830, p = 0.0051); `Fig2_module_scores_donor_UCell.*` (UCell counterpart: 0.407 vs 0.706, Welch p = 0.0094, d = 7.43).

**Fixed?** Yes — corrected figures + audit table + two consistent text/legend options provided (manuscript-ready file).
**Impact on conclusions.** None on the result itself (both algorithms agree in direction; doublet sensitivity shown for both); removes the algorithm/level mismatch between figure and text.
**Remaining limitations.** The final choice of which version to embed in Figure 2 rests with the authors; legend text for both options is supplied.
**Evidence.** `tables/R2_figure2_scoring_audit.csv`; `logs/R2_task5_figure2_scoring_audit.log`; script `R2_task5_figure2_scoring_audit.R`.

---

## Consolidated verdict

- Annotation contradiction: **resolved and corrected** (stage structure documented; claim fixed).
- Code repository: **reconstructed with executed validation**; local git tag prepared; nothing published online.
- External negative control: **empirical-specificity conclusions stand**; direction-pattern and case-control vs empirical distinctions now stated correctly.
- Meta-analysis: **reproduced exactly**; pooled random-effects test for Fibrocartilage is not significant (p = 0.107) — descriptive synthesis supplied and inferential claims must be softened.
- Figure 2: **algorithm/level mismatch documented and consistent figures supplied**.
- The paper's primary donor-level finding (Fibrocartilage program higher in ONFH; unchanged after doublet removal) is intact and now has a complete, honest evidence chain.
