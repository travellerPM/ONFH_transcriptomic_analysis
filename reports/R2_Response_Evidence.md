# R2 Response Evidence — Reviewer-2 revision support material

**Note on use.** Reviewer-2 comments are listed below by topic as handled in this revision; the authors should paste the verbatim comment text before submission. No manuscript line numbers or "location of changes" are invented here. All numbers are taken from the audited result files listed per item. Draft responses are suggestions for the authors to edit.

---

## Item 1 — Donor-level statistics (cells are not independent replicates)

- **Analysis performed:** donor-level aggregation of program scores; Welch t-test and Cohen's d at n = 3 vs 3; donor-level figure.
- **Exact results:** Fibrocartilage donor means (AddModuleScore MS_2): Control 0.4275/0.7067/0.5333 (group mean 0.556 ± 0.141) vs ONFH 1.5212/1.9831/1.9843 (1.830 ± 0.267); Welch t = 7.306, p = 0.0051, Cohen's d = 5.97. UCell: 0.407 ± 0.011 vs 0.706 ± 0.056, p = 0.0094, d = 7.43.
- **Supporting files:** `tables/donor_module_scores.csv`, `tables/R2_donor_scores_group_stats.csv`, `tables/R2_donor_scores_all_programs_stats.csv`, `figures/R2_Fig_donor_scores_corrected.svg/.png` (label-corrected), `figures/Fig2_module_scores_donor_AddModuleScore.svg/.png`, `figures/Fig2_module_scores_donor_UCell.svg/.png`.
- **Suggested response (draft):** "All program-score comparisons are performed at the donor level, with donors as the biological unit; single-cell values are shown only as distributions. We have corrected the axis label of the donor-level figure (the data are AddModuleScore scores) and provide a donor-level figure panel that matches the quoted statistics exactly."
- **Remaining limitations:** n = 3 per group; minimum two-sided exact permutation p = 0.10; parametric tests assume normality.

## Item 2 — Cell composition and canonical marker confirmation

- **Analysis performed:** per-sample composition (total, cartilage-lineage, fibrocartilage-candidate, fibrocartilage-dominant); per-sample SOX9/ACAN/COL2A1 detection for each ONFH sample vs other cartilage-lineage cells (descriptive only).
- **Exact results:** composition A-01 5,892/5,757/135/1; A-03 15,838/14,880/958/13; A-07 9,458/9,160/298/2; A-02 9,225/77/9,148/2,678; A-06 10,023/34/9,989/2,988; A-08 6,874/7/6,867/1,467. ACAN detection higher in Fibrocartilage-dominant cells in all three ONFH samples (76.2 vs 35.5%; 72.2 vs 32.3%; 81.5 vs 57.1%); comparator groups are small (76/31/7 cells) and flagged as unstable.
- **Supporting files:** `tables/R2_cell_composition_main_table.csv`, `tables/R2_canonical_markers_per_ONFH_sample.csv`, `tables/R2_canonical_markers_per_sample_all.csv`, `figures/R2_Fig_canonical_markers_per_sample.svg/.png`, `R2_task2_task3_completion_report.md`.
- **Suggested response (draft):** "We now report per-sample composition and per-sample canonical-marker detection. Chondrocyte markers remain detectable in Fibrocartilage-dominant cells in every ONFH sample (ACAN consistently higher than in other cartilage-lineage cells), supporting a fibrocartilage-like state of cartilage-lineage cells. All per-sample statistics are descriptive; no significance tests were applied to the small comparator groups."
- **Remaining limitations:** comparator groups in ONFH samples contain 76/31/7 cells; Control Fibrocartilage-dominant groups contain 1–13 cells.

## Item 3 — Doublet confounding

- **Analysis performed:** scDblFinder (per-donor); doublet fractions in Fibrocartilage-dominant (A) vs other cartilage-lineage (B) cells with the reviewer-requested comparator definition; donor-level paired Wilcoxon; condition-balanced sensitivity (A vs Fibrocartilage-candidate C; A vs B within ONFH); donor-level scores before/after doublet removal (AddModuleScore primary, UCell secondary) with Welch tests, Cohen's d, exact permutation (all 20 splits) and Monte-Carlo permutation (100,000 draws, seed 42).
- **Exact results:** A 8.09% (578/7,149) vs B 11.08% (3,313/29,899), Fisher p = 2.7e-14, OR = 0.706 (descriptive; A is 99.8% ONFH and B 99.6% Control cells); A vs C 8.09% vs 8.91% (p = 0.034, OR = 0.899); donor-level paired Wilcoxon p = 0.84; within ONFH donors the direction reverses (A 8.08% vs B 0.88%, n = 114). Scores before vs after doublet removal: MS_2 0.556 ± 0.141 vs 1.830 ± 0.267 (p = 0.0051, d = 5.97) → 0.555 ± 0.139 vs 1.822 ± 0.269 (p = 0.0054, d = 5.92); UCell p = 0.0094 → 0.0085. Exact permutation one-sided 0.05 / two-sided 0.10; Monte-Carlo ≈ 0.10.
- **Supporting files:** `tables/R2_dblsens_cartilage_lineage_rates.csv`, `..._per_donor.csv`, `..._score_sensitivity.csv`, `..._fraction_tests.csv`, `..._donor_scores.csv`, `tables/R2_dblsens_condition_sensitivity.csv`, `figures/R2_Fig_dblsens_cartilage_lineage.svg/.png`, `R2_cell_identity_audit.md`.
- **Suggested response (draft):** "We rebuilt the doublet comparator exactly as requested (Fibrocartilage-dominant vs other cartilage-lineage cells per the original lineage_call definition). Predicted doublets are not enriched in Fibrocartilage-dominant cells, the cell-level difference is confounded by disease composition and is not significant at the donor level, and removing predicted doublets leaves the donor-level Fibrocartilage difference essentially unchanged. We therefore consider the program difference robust to doublet effects."
- **Remaining limitations:** scDblFinder calls are predictions; Control Fibrocartilage-dominant groups contain 1–13 cells; minimum two-sided permutation p = 0.10 at n = 3 vs 3.

## Item 4 — Negative control / signature specificity in external cohorts

- **Analysis performed:** ssGSEA of the frozen Fibrocartilage signature vs 1,000 expression-matched random sets and EMT/ECM/fibroblast comparators in GSE74089 and GSE345158; statistical review of the observed-vs-null definitions with symmetric correction; standardized effect sizes.
- **Exact results:** disease effect vs random null: GSE74089 0.400, empirical p = 0.001, z = 5.37; GSE345158 0.122, p = 0.002, z = 4.64 (definition already symmetric; recomputation identical). Directional concordance: observed 0.957 (22/23) vs symmetric null 0.577, corrected p = 0.0010 (GSE74089); 0.739 (17/23) vs 0.387, corrected p = 0.0030 (GSE345158); direction-permutation null p = 1.0 / 0.326 — the frozen direction pattern itself is not informative externally (DCN down-direction not reproduced). Case-control t-tests: GSE74089 p = 0.0025; GSE345158 p = 0.253 (not significant). Cohen's d: GSE74089 Fibrocartilage 6.58 vs EMT 9.47 / ECM 14.39 / Fibroblast 2.74; GSE345158 0.90 vs −0.22 / 0.13 / 0.82.
- **Supporting files:** `tables/R2_extspec_random_results.csv`, `R2_extspec_random_null_raw.csv`, `R2_extspec_comparator_results.csv`, `R2_extspec_stats_review.csv`, `R2_extspec_corrected_concordance.csv`, `figures/R2_Fig_extspec.svg/.png`.
- **Suggested response (draft):** "We reviewed the concordance statistic and replaced the up-fraction null with a symmetric null that assigns the frozen direction pattern to matched random sets; the empirical-specificity conclusions are unchanged (p = 0.001 / 0.003). We now report three distinct quantities — empirical significance vs random sets, case-control significance, and relative specificity vs generic programs — and standardized effect sizes because raw ssGSEA effects are not comparable across signature sizes. We explicitly note that the concordance reflects coordinated up-regulation of signature genes in disease and that GSE345158 alone is not significant at case-control level."
- **Remaining limitations:** small external cohorts (4v4, 5v4); comparator-specific matched nulls not generated (t/d/z used instead); ssGSEA scores not variance-homogenized across platforms.

## Item 5 — Independent re-run of GSE74089

- **Analysis performed:** independent re-derivation of probe annotation and limma DE from the raw series matrix.
- **Exact results:** Spearman rho(logFC) = 1.0 vs the original analysis; 8,952/8,952 significant genes reproduced.
- **Supporting files:** `tables/R2_gse74089_rerun_verification.csv`, script `R2_task5_gse74089_rerun.R`.
- **Suggested response (draft):** "An independent re-run from the raw files reproduced the original GSE74089 differential-expression results exactly (rho = 1.0)."
- **Remaining limitations:** none material.

## Item 6 — Code and parameter availability

- **Analysis performed:** exhaustive search for the lost Phase-1 scripts; reconstruction of eight equivalent executable scripts (all headed RECONSTRUCTED, UNKNOWN parameters flagged); executed consistency checks against persisted results; README update; local version control with annotated tag.
- **Exact results:** V1 AddModuleScore recompute reproduces stored scores exactly (per-cell Spearman = 1.000, all 5 programs; donor-mean diff = 0.0000); V2 edgeR pseudobulk recompute rho(logFC) = 0.9999 (18,450 genes; FDR<0.05 Jaccard 0.826, consistent with version drift); V3 res 0.4 = 13 clusters, no batch integration; V4 QC reconciliation exact; V5 AUCell columns unpopulated as documented. GSE74089 re-run rho = 1.0; ssGSEA reproduction rho = 0.982 / 0.787.
- **Supporting files:** `code_repository/` (README, parameters.md, sessionInfo, all scripts, `scripts/reconstructed_phase1/`), `tables/RECONSTRUCTED_validation_results.csv`, `logs/RECONSTRUCTED_validation_check.log`, `RELEASE_NOTES_v2.0-revision-r2.md`.
- **Suggested response (draft):** "The complete analysis code is provided. Original session-local copies of seven Phase-1 scripts were lost in a storage reset; equivalent scripts were rebuilt from the persisted objects and documented parameters, are explicitly marked RECONSTRUCTED, and were validated against the saved results (module-score reproduction exact; pseudobulk logFC rho = 0.9999). Steps that cannot be re-executed are marked unverified. The package is version-controlled (tag v2.0-revision-r2); the tag is local and has not been pushed to a public repository."
- **Remaining limitations:** raw-QC thresholds, integration-free clustering details, and the phase-0.5 marker-gene lists could not be recovered and are marked UNKNOWN; robustness re-execution not performed (persisted tables authoritative).

## Item 7 — Cell annotation provenance (raised by the integrity audit)

- **Analysis performed:** barcode-level tracking across persisted objects; stage-by-stage count reconciliation; per-cell cross-tabulation; group-partition audit; phase-0 panel-logic reconstruction; condition-composition sensitivity analysis.
- **Exact results:** 57,476 QC'd cells = 57,310 phase-0.5 Chondrocyte + 166 Macrophage-cluster cells (removed; per-donor reconciliation 6/6 exact); A+B+C = 57,310 complete partition with zero overlap; A is 99.78% ONFH cells; 97.0% of A cells co-express chondrocyte and fibroblast-panel markers; A-vs-C doublet comparison 8.09% vs 8.91% (donor-level p = 0.59).
- **Supporting files:** `R2_cell_identity_audit.md`, `tables/R2_cell_annotation_provenance.csv`, `R2_lineage_crosswalk.csv`, `R2_cell_identity_crosstab.csv`, `R2_group_composition_audit.csv`, `R2_phase0_panel_score_by_group.csv`, `R2_marker_coexpression_by_group.csv`, `R2_dblsens_condition_sensitivity.csv`.
- **Suggested response (draft):** "We clarified the annotation provenance: the analysed cells are cartilage-lineage cells retained after marker-based lineage re-annotation of the integrated object; 166 cells forming a macrophage-marker cluster were excluded, and the count reconciliation is exact at every stage. We corrected an imprecise statement in an internal report that had described the retained cells as 'phase-0 chondrocyte calls'. We further added a sensitivity analysis showing that the cell-level doublet comparison between Fibrocartilage-dominant and other cartilage-lineage cells is confounded by disease composition, and we now cite the donor-level and condition-balanced results instead."
- **Remaining limitations:** per-cell phase-0 labels and the 166 removed barcodes are not recoverable from persisted artifacts.

## Item 8 — Pooled external evidence (meta-analysis)

- **Analysis performed:** exact recomputation of the DerSimonian-Laird synthesis; Hedges-g SMD sensitivity; fixed-effect comparison; descriptive synthesis.
- **Exact results:** recomputation matches the published table exactly (max diff = 0). Fibrocartilage: DL combined 0.197 [−0.042, 0.436], p = 0.107, I² = 93.1%, τ² = 0.028; FE 0.247, p = 9.3e-18 (invalid under heterogeneity); Hedges g 5.72 (GSE74089) vs 0.836 (GSE345158); SMD-DL 3.00 [−1.76, 7.75], p = 0.216. Inflammatory direction-discordant (−0.315 vs +0.050).
- **Supporting files:** `tables/R1_meta_analysis_results.csv` (original), `tables/R1_meta_audit_results.csv` (audit), `logs/R2_task4_meta_analysis_audit.log`, script `R2_task4_meta_analysis_audit.R`.
- **Suggested response (draft):** "We audited the random-effects synthesis. The model, effect metric, and standard errors are reproduced exactly. Because pooling two studies leaves one degree of freedom for heterogeneity, we now present the pooled estimate descriptively alongside the per-cohort effects and direction concordance, report a standardized (Hedges g) sensitivity analysis, and do not interpret I²/τ² substantively. The pooled random-effects test for the Fibrocartilage signature is not significant (p = 0.107); external support rests on the direction-concordant, individually significant GSE74089 effect and the concordant direction in GSE345158."
- **Remaining limitations:** two cohorts; SMD sensitivity only for signatures with saved per-sample scores; small per-cohort sample sizes.
