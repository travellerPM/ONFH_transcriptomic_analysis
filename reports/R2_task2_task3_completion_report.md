# R2 Task 2 + Task 3 Completion Report
## Doublet sensitivity with cartilage-lineage comparator & per-sample canonical marker analysis (GSE311084)

Date: 2026-10-15
Scope: Only the two remaining Reviewer-2 analysis tasks were executed. The external negative control, the random-effects meta-analysis, the Abstract, Figure 2, and the README were **not** modified or re-run. No Response to Reviewers was drafted. All previously generated results (including the first-version `R2_dblsens_*` files) are preserved unchanged.

---

## 1. Environment and reproducibility

| Item | Value |
|---|---|
| R version | 4.4.2 (2024-10-31) |
| Key packages | Seurat 5.5.1, SeuratObject 5.4.0, ggplot2 4.0.3, patchwork 1.3.2, svglite 2.2.2, dplyr 1.2.1 |
| Random seed | 42 (set at the top of both scripts; used for the Monte-Carlo permutation) |
| scDblFinder classifications | Pre-computed `scDblFinder_class` columns in `GSE311084_scored_doublet.rds` (from the earlier scDblFinder run); scDblFinder was **not** re-run |
| UCell scores | Pre-computed `FibrocartilageUC_` columns in the RDS metadata; UCell was **not** re-run |
| Console logs | `/mnt/results/phase1/logs/R2_task3_doublet_cartilage_lineage_comparator.log`, `/mnt/results/phase1/logs/R2_task2_per_sample_canonical_markers.log` |
| sessionInfo | `R2_dblsens_cartilage_lineage_sessionInfo.txt`, `R2_canonical_markers_per_sample_sessionInfo.txt` (in tables/) |

Both scripts exited with status 0 and full console logs were captured.

---

## 2. Comparator definitions (shared by Task 1 and Task 2)

**A — Fibrocartilage-dominant**
```r
dominant_program == "Fibrocartilage"
```

**B — Other cartilage-lineage cells (primary comparator)**
```r
lineage_call == "cartilage_lineage" & dominant_program != "Fibrocartilage"
```

**Provenance of `lineage_call` (recovered from the phase-0.5 script and `GSE311084_lineage.rds`):** `lineage_call` is a *cluster-level* call made from supportive marker counts across 9 lineages. Clusters with 5 chondrocyte markers and 0 fibroblast/stromal markers were called `cartilage_lineage`; clusters with **both** chondrocyte and fibroblast/stromal markers were called `cartilage_fibrocartilage_candidate`. All 57,310 GSE311084 cells are Chondrocyte calls from the **phase-0.5 integrated re-annotation** — there are no macrophage, undetermined-lineage, or non-cartilage-lineage cells in the retained set, so the primary comparator B is not contaminated by non-cartilage cells.

> **Correction (2026-10-15, R2 integrity audit):** an earlier version of this report stated that all cells are "phase-0 Chondrocyte calls". That was imprecise: phase-0 *per-sample preliminary* annotations were heterogeneous (Fibroblast 38,585 / Chondrocyte 18,765 / Macrophage 126 of 57,476 QC'd cells), and the 166 cells in the phase-0.5 Macrophage cluster were removed before the final 57,310-cell object. The full provenance, barcode-level tracking, and a condition-composition sensitivity analysis for the A-vs-B doublet comparison are documented in `R2_cell_identity_audit.md`, `tables/R2_cell_annotation_provenance.csv`, and `tables/R2_dblsens_condition_sensitivity.csv`.

**Note on `dominant_program == "Undetermined"` (8,894 cells within B):** this label denotes the *absence of a dominant remodeling program* (program-level annotation), **not** lineage uncertainty. These cells are retained in the primary comparator B. As a documented sensitivity, **B2 — B without Undetermined** (`lineage_call == "cartilage_lineage" & !dominant_program %in% c("Fibrocartilage","Undetermined")`, n = 21,005) was analysed in parallel.

**Context groups reported for transparency (not primary comparators):**
- C — Other fibrocartilage-candidate cells: `lineage_call == "cartilage_fibrocartilage_candidate" & dominant_program != "Fibrocartilage"` (n = 20,262)
- REF — previous first-version definition "all non-fibrocartilage cells" (n = 50,161), included only as a reference row so old and new definitions can be compared directly.

**Missing-data policy:** if a donor had zero cells in a group, the doublet percentage / marker statistics would be reported as `NA` together with the actual denominator (n_cells = 0). No imputation was performed. In practice no donor×group cell was empty; the smallest strata are flagged as unstable instead (see below).

---

## 3. Task 1 — Doublet sensitivity with the cartilage-lineage comparator

### 3.1 Group-level doublet rates (`R2_dblsens_cartilage_lineage_rates.csv`)

| Group | n cells | n doublets | doublet % |
|---|---|---|---|
| A Fibrocartilage-dominant | 7,149 | 578 | 8.09 |
| B Other cartilage-lineage | 29,899 | 3,313 | 11.08 |
| B2 B without Undetermined | 21,005 | 2,580 | 12.28 |
| C Other fibrocartilage-candidate | 20,262 | 1,806 | 8.91 |
| REF previous "Other cells" (reference only) | 50,161 | 5,119 | 10.21 |

### 3.2 Per-donor rates (`R2_dblsens_cartilage_lineage_per_donor.csv`)

| Donor | A: n / doublet % | B: n / doublet % |
|---|---|---|
| Ctrl-1 | 1 / 0.00 | 5,756 / 6.18 |
| Ctrl-2 | 13 / 15.38 | 14,870 / 13.32 |
| Ctrl-3 | 2 / 0.00 | 9,159 / 10.66 |
| ONFH-1 | 2,678 / 9.41 | 76 / 1.32 |
| ONFH-2 | 2,988 / 8.43 | 31 / 0.00 |
| ONFH-3 | 1,467 / 4.91 | 7 / 0.00 |

### 3.3 Fraction tests (`R2_dblsens_cartilage_lineage_fraction_tests.csv`)

| Test | Result |
|---|---|
| Fisher exact, A vs B (cell-level, **descriptive**) | p = 2.67e-14, OR = 0.706 |
| Fisher exact, A vs B2 (sensitivity) | p = 1.42e-23, OR = 0.628 |
| Donor-level paired Wilcoxon, A vs B (n = 6 pairs) | p = 0.844 (n.s.) |

**Interpretation.** Predicted doublets are modestly **depleted**, not enriched, in Fibrocartilage-dominant cells (OR < 1 in both comparisons). The donor-level paired test is not significant, and it is strongly limited here: the Control A-groups contain only 1, 13 and 2 cells (Ctrl-1/2/3), and Ctrl-2 shows the opposite direction (15.38% vs 13.32%). The cell-level Fisher p-value is descriptive only — cells are not independent replicates; the donor is the biological unit.

### 3.4 Score sensitivity before vs after doublet removal (`R2_dblsens_cartilage_lineage_score_sensitivity.csv`, donor scores in `R2_dblsens_cartilage_lineage_donor_scores.csv`)

Donor-level fibrocartilage scores were recomputed on all cells and on singlets only (n_all / n_singlet per donor recorded). Group statistics (n = 3 vs 3 donors):

| Score | Analysis | Control mean ± SD | ONFH mean ± SD | Welch t | Welch p | Cohen's d |
|---|---|---|---|---|---|---|
| AddModuleScore `MS_2` (primary) | before | 0.556 ± 0.141 | 1.830 ± 0.267 | 7.306 | 0.0051 | 5.97 |
| AddModuleScore `MS_2` (primary) | after | 0.555 ± 0.139 | 1.822 ± 0.269 | 7.252 | 0.0054 | 5.92 |
| UCell (secondary) | before | 0.407 ± 0.011 | 0.706 ± 0.056 | 9.105 | 0.0094 | 7.43 |
| UCell (secondary) | after | 0.404 ± 0.010 | 0.700 ± 0.053 | 9.620 | 0.0085 | 7.85 |

Permutation tests on the `MS_2` donor-level difference (Control vs ONFH):
- Exact permutation (all C(6,3) = 20 donor-label splits enumerated): one-sided p = 0.05, two-sided p = 0.10 — the minimum attainable two-sided p at n = 3 vs 3.
- Monte-Carlo permutation (100,000 draws, seed 42): two-sided p = 0.101 (before) / 0.099 (after).

**Interpretation.** Removing predicted doublets leaves the donor-level fibrocartilage score difference essentially unchanged (direction unchanged; d 5.97 → 5.92; Welch p 0.0051 → 0.0054). The result is not driven by doublets. With 3 donors per group, the exact two-sided permutation p is bounded at 0.10; the Welch p-values are reported for continuity with the manuscript but rest on the normality assumption.

---

## 4. Task 2 — Per-sample canonical marker analysis (SOX9 / ACAN / COL2A1)

### 4.1 Per-sample composition (`R2_cell_composition_main_table.csv`)

| Sample | Donor | Condition | Total | cartilage_lineage | fibrocartilage_candidate | fibrocartilage_dominant | % FC-dominant |
|---|---|---|---|---|---|---|---|
| A-01 | Ctrl-1 | Control | 5,892 | 5,757 | 135 | 1 | 0.02 |
| A-03 | Ctrl-2 | Control | 15,838 | 14,880 | 958 | 13 | 0.08 |
| A-07 | Ctrl-3 | Control | 9,458 | 9,160 | 298 | 2 | 0.02 |
| A-02 | ONFH-1 | ONFH | 9,225 | 77 | 9,148 | 2,678 | 29.03 |
| A-06 | ONFH-2 | ONFH | 10,023 | 34 | 9,989 | 2,988 | 29.81 |
| A-08 | ONFH-3 | ONFH | 6,874 | 7 | 6,867 | 1,467 | 21.34 |

Sample IDs in the metadata match the expected mapping (A-02/ONFH-1, A-06/ONFH-2, A-08/ONFH-3); metadata was used as truth.

### 4.2 ONFH per-sample marker detection (`R2_canonical_markers_per_ONFH_sample.csv`; all 6 samples in `R2_canonical_markers_per_sample_all.csv`, 36 rows)

Percent expressing (count > 0 on log-normalized data), Fibrocartilage-dominant (A) vs Other cartilage-lineage (B):

| Sample | Group (n) | SOX9 | ACAN | COL2A1 |
|---|---|---|---|---|
| ONFH-1 | A (2,678) | 63.5% | 76.2% | 45.2% |
| ONFH-1 | B (76) | 61.8% | 35.5% | 36.8% |
| ONFH-2 | A (2,988) | 46.2% | 72.2% | 55.4% |
| ONFH-2 | B (31) | 67.7% | 32.3% | 45.2% |
| ONFH-3 | A (1,467) | 50.9% | 81.5% | 99.2% |
| ONFH-3 | B (7) | 42.9% | 57.1% | 100% |

**Conventions.** `mean_log_normalized_expression` and `median_log_normalized_expression` are computed over **all cells in the group, zeros included** (Seurat DotPlot convention), so they jointly reflect detection frequency and expression level. Detection = n_expressing / n_cells × 100.

**Key observations (descriptive only).**
- ACAN detection is higher in Fibrocartilage-dominant cells in **all three** ONFH samples (76.2 vs 35.5; 72.2 vs 32.3; 81.5 vs 57.1).
- SOX9 and COL2A1 detection are variable across samples and do not show a consistent direction.
- Canonical chondrocyte markers remain detectable in Fibrocartilage-dominant cells in every ONFH sample, supporting chondrocyte identity of these cells at the per-sample level.

**Small-sample flags (recorded in the `note` column).** The B comparator is very small in ONFH samples (ONFH-1 n = 76, ONFH-2 n = 31, ONFH-3 n = 7); ONFH-3 B rows are explicitly flagged "small group (n=7); descriptive value unstable". Control A-groups are n = 1 / 13 / 2 and are flagged the same way in the all-samples table. No significance tests were run per sample; all statistics are descriptive. No group was empty, so no NA cells were required, but the NA + denominator policy is implemented in the script.

---

## 5. Output files

**Scripts**
- `/workspace/phase1/scripts/R2_task3_doublet_cartilage_lineage_comparator.R` (Task 1)
- `/workspace/phase1/scripts/R2_task2_per_sample_canonical_markers.R` (Task 2)
- Copies added to the code repository: `code_repository/scripts/R2_task3_doublet_cartilage_lineage_comparator.R`, `code_repository/scripts/R2_task2_per_sample_canonical_markers.R` (README not modified)

**Tables** (`/mnt/results/phase1/tables/`)
- `R2_dblsens_cartilage_lineage_rates.csv`
- `R2_dblsens_cartilage_lineage_per_donor.csv`
- `R2_dblsens_cartilage_lineage_score_sensitivity.csv`
- `R2_dblsens_cartilage_lineage_fraction_tests.csv`
- `R2_dblsens_cartilage_lineage_donor_scores.csv` (per-donor MS_2/UCell before vs after; supplementary)
- `R2_dblsens_cartilage_lineage_sessionInfo.txt`
- `R2_canonical_markers_per_ONFH_sample.csv` (18 rows)
- `R2_canonical_markers_per_sample_all.csv` (36 rows)
- `R2_cell_composition_main_table.csv`
- `R2_canonical_markers_per_sample_sessionInfo.txt`

**Figures** (`/mnt/results/phase1/figures/`)
- `R2_Fig_dblsens_cartilage_lineage.svg` / `.png`
- `R2_Fig_canonical_markers_per_sample.svg` / `.png`

**Logs** (`/mnt/results/phase1/logs/`)
- `R2_task3_doublet_cartilage_lineage_comparator.log`
- `R2_task2_per_sample_canonical_markers.log`

---

## 6. Preservation and scope statements

- The first-version doublet sensitivity results (`R2_dblsens_sensitivity_summary.csv`, `R2_dblsens_*` files, `R2_Fig_dblsens*`) and all other existing results were **not** deleted or overwritten; the new outputs use distinct `*_cartilage_lineage*` / `*_per_sample*` filenames.
- The external negative control and the random-effects meta-analysis were not re-run; the Abstract, Figure 2, and README were not modified; no Response to Reviewers was drafted.
- All analyses here are **robustness/descriptive assessments of the existing GSE311084 dataset**. No new biological-validation claims are made: the doublet analysis addresses whether the donor-level score difference is confounded by predicted doublets; the per-sample marker analysis documents chondrocyte-marker retention per sample. Cell-level p-values are descriptive; the donor is the biological unit throughout.
- Known limitations: (i) Control A-groups contain 1–13 cells, so donor-level doublet comparisons in Control donors are unstable; (ii) the ONFH B comparator is small (76/31/7 cells); (iii) with n = 3 vs 3 donors the minimum two-sided exact permutation p is 0.10; (iv) scDblFinder classifications are predictions, not experimental ground truth.
