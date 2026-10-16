# R2 Cell Identity Audit — GSE311084 Annotation Provenance

**Date:** 2026-10-15
**Trigger:** Contradiction between (i) `R2_cell_composition_per_sample.csv`, which records phase-0 per-sample preliminary annotations (Fibroblast 38,585 / Chondrocyte 18,765 / Macrophage 126 = 57,476 QC'd cells), and (ii) a comment in `R2_task3_doublet_cartilage_lineage_comparator.R` claiming "All cells in GSE311084 are phase-0 Chondrocyte calls" for the 57,310-cell final analysis set.
**Verdict up front:** the comment was **wrong as stated** and has been corrected. The two descriptions refer to **different annotation stages**. Full provenance below.

---

## 1. Verification method

Artifacts used (all persisted; no re-analysis of raw data was required):

| Stage | Artifact | Content |
|---|---|---|
| Phase 0 (per-sample QC + preliminary annotation) | `onfh_data/phase0_celltype_results.json` | Per-sample: raw/QC'd counts, 15 clusters/sample, cluster-level calls from 13-lineage marker panels, per-cluster mean marker scores |
| Phase 0.5 (integrated re-annotation) | `phase0_5/GSE311084_major_calls.rds`, `GSE311084_lineage.rds`, `GSE311084_ct_counts.rds`, `GSE311084_ct_donor.rds`, `GSE311084_qc_summary.rds` | Cluster-level lineage calls with rationale; aggregate cell-type counts |
| Final analysis objects | `phase1/GSE311084_reclustered.rds`, `GSE311084_scored.rds`, `GSE311084_scored_doublet.rds` | Per-cell metadata: `cell_type`, `lineage_call`, `dominant_program`, `donor_id`, `sample_id`, `condition` |
| Revision table | `tables/R2_cell_composition_per_sample.csv` | Source of the 38,585/18,765/126 figures |

Barcode-level operations were performed in `R2_task1_cell_annotation_provenance.R` (log: `logs/R2_task1_cell_annotation_provenance.log`).

## 2. What the three stages are — and what "Fibroblast/Chondrocyte/Macrophage" means

1. **Phase 0 — per-sample preliminary annotation (NOT the final cell identity).** Each of the 6 samples was QC'd and clustered independently (15 clusters per sample); clusters were called against 13 lineage marker panels (neutrophil, monocyte, macrophage, T cell, NK cell, B cell, endothelial, **fibroblast [COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRA]**, **chondrocyte [COL2A1, ACAN, SOX9, COL9A1/2, COL11A1, HAPLN1]**, osteoblast, osteoclast, pericyte, MSC). Result: Fibroblast 38,585 / Chondrocyte 18,765 / Macrophage 126 of 57,476 QC'd cells. These are **preliminary, per-sample, cluster-level calls**, superseded in phase 0.5.
2. **Phase 0.5 — integrated re-annotation (the lineage framework used by the paper).** Samples were integrated and re-clustered (11 clusters at the lineage-calling step); each cluster was called from supportive-marker counts across 9 lineages (Chondrocyte, Fibroblast_stromal, Endothelial, Monocyte, Macrophage, T cell, Mast cell, Osteoclast, Neutrophil). Persisted rationale strings (`GSE311084_major_calls.rds`): clusters 0–9 → Chondrocyte, with `cartilage_lineage` (cartilage markers present, fibroblast/stromal = 0) vs `cartilage_fibrocartilage_candidate` (both cartilage and fibroblast/stromal markers present) recorded in `lineage_call`; cluster 10 → **Macrophage** ("Competitive: Macrophage has 4 supportive markers vs Chondrocyte 3"), 166 cells (ONFH-1: 7, ONFH-2: 159).
3. **Final analysis set.** The 166 phase-0.5 Macrophage-cluster cells were removed; the remaining **57,310 cells all carry `cell_type = Chondrocyte` from the phase-0.5 re-annotation**, plus per-cell `lineage_call` (cartilage_lineage 29,915; cartilage_fibrocartilage_candidate 27,395) and `dominant_program`.

**Correction issued:** the correct statement is "all 57,310 retained cells are Chondrocyte calls from the **phase-0.5 integrated re-annotation**", not "phase-0 Chondrocyte calls". The comment in `R2_task3_doublet_cartilage_lineage_comparator.R` and the corresponding sentence in `R2_task2_task3_completion_report.md` have been corrected in place with dated correction notes. The terminology was not changed to mask the contradiction; the underlying stage structure is documented here.

## 3. Reconciliation: 57,476 → 57,310

| Sample | Donor | Phase-0 QC'd | Phase-0.5 Chondrocyte | Phase-0.5 Macrophage | Final | Removed | Removed = Macrophage? |
|---|---|---|---|---|---|---|---|
| A-01 | Ctrl-1 | 5,892 | 5,892 | 0 | 5,892 | 0 | TRUE |
| A-02 | ONFH-1 | 9,232 | 9,225 | 7 | 9,225 | 7 | TRUE |
| A-03 | Ctrl-2 | 15,838 | 15,838 | 0 | 15,838 | 0 | TRUE |
| A-06 | ONFH-2 | 10,182 | 10,023 | 159 | 10,023 | 159 | TRUE |
| A-07 | Ctrl-3 | 9,458 | 9,458 | 0 | 9,458 | 0 | TRUE |
| A-08 | ONFH-3 | 6,874 | 6,874 | 0 | 6,874 | 0 | TRUE |
| **Total** | | **57,476** | **57,310** | **166** | **57,310** | **166** | **6/6 exact** |

The 166-cell difference is fully accounted for by the phase-0.5 Macrophage-cluster removal; nothing else was dropped between phase-0 QC and the final object. (`tables/R2_cell_annotation_provenance.csv`)

## 4. Barcode-level tracking — what could and could not be verified

**Verified 1:1 at barcode level:**
- `reclustered` = `scored` = `scored_doublet`: 57,310 barcodes each; pairwise set differences = 0; no duplicated barcodes.
- 10x barcode suffix (1–6) maps 1:1 to sample_id (A-01…A-08); `orig.ident` constant; donor/sample/condition mapping matches the phase-0 JSON exactly (A-01→Ctrl-1 … A-08→ONFH-3).

**Not recoverable from persisted artifacts (stated honestly, not inferred from aggregates):**
- **Per-cell phase-0 labels.** Phase-0 was per-sample clustering; only cluster-level aggregates persist (JSON + the revision CSV). The per-sample phase-0 objects were lost in the workspace reset, so a cell-by-cell phase-0 → phase-0.5 → final barcode trace is **not possible**; any claim of the form "cell X was originally a fibroblast" cannot be verified and is not made.
- **The 166 removed Macrophage-cluster barcodes.** The phase-0.5 intermediate object (11-cluster lineage-calling clustering) was not persisted; `reclustered.rds` already contains only the 57,310 retained cells. Their removal is verified by exact per-donor count reconciliation (Section 3), not by barcode setdiff.

**Rule-based reconstruction of the phase-0 scoring logic (per cell, descriptive):** applying the phase-0 fibroblast vs chondrocyte panel means to the final object (`tables/R2_phase0_panel_score_by_group.csv`): 100.0% of group A cells, 94.5% of B, 68.1% of C have higher mean expression on the phase-0 fibroblast panel than on the chondrocyte panel. This shows the phase-0 "fibroblast" majority was a **panel artifact**: DCN, LUM, COL1A1, COL1A2, COL3A1 are fibrocartilage ECM genes highly expressed by cartilage-lineage cells in this tissue, and the phase-0 panel had no fibrocartilage category. It is a reconstruction of the rule, not a recovery of the original per-cell calls.

## 5. Cross-tabulation and group audit

- Full per-cell cross-tab (`cell_type × lineage_call × dominant_program × condition × donor_id`): `tables/R2_cell_identity_crosstab.csv` (67 non-zero combinations). A 5-way table including phase-0 calls is not producible (Section 4); the cluster-level crosswalk across all three stages is in `tables/R2_lineage_crosswalk.csv` (90 phase-0 sample-clusters, 11 phase-0.5 clusters, 102 final res-0.4 combinations).
- **Group definitions are correct and non-overlapping:** A (Fibrocartilage-dominant, 7,149) + B (other cartilage_lineage, 29,899) + C (other candidate, 20,262) = 57,310 — a complete partition; A∩B = A∩C = B∩C = 0; B2 (21,005) ⊂ B.
- **Composition imbalance (material finding):** A is 99.78% ONFH cells (7,133/7,149) and 99.78% `cartilage_fibrocartilage_candidate`; B is 99.62% Control cells (29,785/29,899) and 100% `cartilage_lineage`. The cell-level A-vs-B doublet comparison is therefore almost fully confounded with disease condition.

## 6. Does group A contain "originally Fibroblast" cells — and what it means biologically

Per-cell phase-0 labels are not recoverable (Section 4), so this cannot be proven cell-by-cell. Two facts bound the answer:
1. By aggregate counts, 67% of all QC'd cells (38,585/57,476) were phase-0 fibroblast calls, including 64–86% of cells in the ONFH samples that contribute 99.8% of group A. It is therefore near-certain that most A-group cells carried a preliminary phase-0 "fibroblast" call.
2. Per-cell marker evidence in the final object: A-group cells co-express both panels — 97.0% express ≥1 chondrocyte-panel gene (SOX9 53.6%, ACAN 75.5%, COL2A1 60.5%) and 100% express ≥1 fibroblast-panel gene; 97.0% express both (`tables/R2_marker_coexpression_by_group.csv`).

**Biological interpretation (unchanged, now properly grounded):** the Fibrocartilage-dominant cells are cartilage-lineage cells (phase-0.5 marker-based calls; canonical chondrocyte markers retained per sample) that co-express fibroblastic/fibrocartilage ECM genes — which is the definition of a fibrocartilage phenotype and of the frozen 24-gene signature (COL1A1/2, COL3A1, FMOD, DCN, LUM + data-driven fibroblastic markers). The phase-0 "fibroblast" calls reflect panel limitations, not a different lineage. The manuscript should describe the annotation as marker-based lineage re-annotation of an integrated cartilage object, and must not describe the 57,310 cells as "phase-0 chondrocytes".

## 7. Minimal sensitivity analysis for the A/B composition imbalance

`R2_task1b_sensitivity_condition_composition.R` → `tables/R2_dblsens_condition_sensitivity.csv` (descriptive; donor = biological unit):

| Analysis | Result |
|---|---|
| S1: A vs C (both ONFH-dominant; condition-balanced program comparison), cell-level Fisher | A 8.09% (7,149) vs C 8.91% (20,262); p = 0.034, OR = 0.899 — far smaller than A vs B (OR = 0.706, p = 2.7e-14) |
| S1: A vs C, donor-level paired Wilcoxon | n = 6 pairs, p = 0.59 (n.s.) |
| S2: A vs B within ONFH donors only | A 8.08% (7,133) vs B 0.88% (114 cells); p = 0.0014, OR = 9.9 — **direction reverses**; B(ONFH) is too small for stable estimates |
| S3: condition effect within strata | cartilage_lineage: ONFH 0.85% vs Control 11.12% (p = 2.5e-05); candidate: ONFH 9.05% vs Control 2.23% (p = 2.3e-24) |

**Conclusion:** the cell-level A-vs-B doublet "depletion" is driven by condition/stratum composition, not by the Fibrocartilage program. Program-matched comparisons (A vs C) show no donor-level doublet difference. This is consistent with, and strengthens, the donor-level paired result (p = 0.844) and the before/after score-sensitivity result (donor-level MS_2 difference essentially unchanged by doublet removal: Welch p 0.0051 → 0.0054, d 5.97 → 5.92). For the manuscript, the doublet robustness argument should rest on the **donor-level** tests and the score-sensitivity analysis; the cell-level Fisher values should be quoted only as descriptive statistics with the composition caveat.

## 8. Impact on the paper's conclusions

- **No impact on the primary result.** The donor-level Fibrocartilage score difference (Control 0.556 ± 0.141 vs ONFH 1.830 ± 0.267; Welch p = 0.0051) is a donor-level analysis and is unaffected by the annotation-stage clarification. Doublet removal does not attenuate it.
- **Wording corrections required (implemented in audit deliverables):** (i) "57,310 chondrocytes" should be "57,310 cells retained after marker-based lineage re-annotation (phase-0.5), all called cartilage-lineage; 166 macrophage-cluster cells removed"; (ii) the doublet comparator description should note the condition-composition imbalance and cite the donor-level and A-vs-C sensitivity results.
- **No change** to frozen signatures, sample sets, dominant-program definitions, or the primary hypothesis.

## 9. Remaining limitations

1. Per-cell phase-0 labels and the 166 removed barcodes are not recoverable; phase-0 → final tracking is aggregate-level only (exact per-donor reconciliation).
2. The phase-0.5 lineage-calling clustering (11 clusters) was an intermediate object; only its cluster-level calls persist. The per-cell `lineage_call`/`cell_type` in the final object are the persisted derivatives of that step.
3. The per-cell phase-0 panel comparison is a rule reconstruction, not the original per-sample cluster calls.
4. scDblFinder calls are predictions; all doublet statistics are descriptive.

## 10. Evidence files

| File | Content |
|---|---|
| `scripts/R2_task1_cell_annotation_provenance.R` (+ code repo copy) | Audit script (EXIT 0) |
| `scripts/R2_task1b_sensitivity_condition_composition.R` (+ code repo copy) | Sensitivity script (EXIT 0) |
| `logs/R2_task1_cell_annotation_provenance.log`, `logs/R2_task1b_sensitivity_condition_composition.log` | Full console logs |
| `tables/R2_cell_annotation_provenance.csv` | Stage-by-stage per-sample counts + reconciliation |
| `tables/R2_lineage_crosswalk.csv` | Cluster-level crosswalk across 3 stages |
| `tables/R2_cell_identity_crosstab.csv` | Per-cell 5-way cross-tab (final metadata) |
| `tables/R2_group_composition_audit.csv` | A/B/C sizes, condition and lineage composition |
| `tables/R2_phase0_panel_score_by_group.csv`, `tables/R2_marker_coexpression_by_group.csv` | Phase-0 panel rule reconstruction; marker co-expression |
| `tables/R2_dblsens_condition_sensitivity.csv` | S1–S3 sensitivity results |
| Corrected: `scripts/R2_task3_doublet_cartilage_lineage_comparator.R`, `R2_task2_task3_completion_report.md` | Dated correction notes |
