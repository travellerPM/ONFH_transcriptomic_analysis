# Analysis Parameters (Complete Record)

All parameters used across Phase 1 (discovery/validation) and Phase 2 (revision). Sources: `FINAL_PHASE1_REPORT.md`, `resolution_selection_summary.md`, executed revision scripts, and persisted output tables.

## 1. Discovery cohort GSE311084 (scRNA-seq)

| Parameter | Value |
|-----------|-------|
| Donors | 6 (3 Control: Ctrl-1/2/3; 3 ONFH: ONFH-1/2/3) |
| Cells (chondrocytes, post-QC) | 57,310 (Ctrl-1=5,892; Ctrl-2=15,838; Ctrl-3=9,458; ONFH-1=9,225; ONFH-2=10,023; ONFH-3=6,874) |
| Clustering resolution | 0.4 (13 clusters; sweep 0.2–0.8; all clusters multi-donor, marker-coherent) |
| Program scoring | AddModuleScore + UCell (AUCell omitted: memory at 57k cells) |
| Pseudobulk DE | edgeR (primary); DESeq2, limma-voom (concordance); donor-level, n=3 vs 3 |
| TF activity | decoupleR wmean + dorothea regulons (confidence A–C, 271 TFs) on pseudobulk |

## 2. Cross-context cohort GSE290411 (scRNA-seq)

| Parameter | Value |
|-----------|-------|
| Donors | 4 SIONFH (SIONFH-1=23, SIONFH-2=2,003, SIONFH-3=86, SIONFH-4=1 cells) |
| Cells used | 2,113 |
| Clustering resolution | 0.3 (6 clusters) |
| Critical caveat | 94.8% of cells from SIONFH-2 → all states labeled donor-specific/low-confidence |

## 3. Frozen signatures (8 signatures, 215 genes)

Defined from GSE311084 BEFORE examining validation data (`phase1_07_freeze_signatures.R`, seed 42):

| Signature | Genes | Basis |
|-----------|-------|-------|
| Homeostatic | 11 curated | SOX9, ACAN, COL2A1, COL9A1/2, COL11A1/2, HAPLN1, CILP, CLEC3A, MATN4 |
| Fibrocartilage | 6 curated | COL1A1, COL1A2, COL3A1, FMOD, DCN, LUM |
| Inflammatory | 26 curated | IL6, CXCL1-3, CXCL8, CCL2/5, ISG15, IFIT1-3, IFI6, IFI44, NFKB1/2, STAT1, IRF1/7, GBP1/4/5, CXCL10, OAS1/2, MX1/2 |
| Pre_hypertrophic | 9 curated | CILP2, IL11, MMP3, CXCL3, IHH, PTHLH, BMP2, BMP4, WNT5A |
| Hypertrophic_osteogenic | 7 present | COL10A1, MMP13, RUNX2, SPP1, MMP11, ALPL, IBSP (BGLAP absent from panel) |
| Conserved_ONFH_remodeling | composite | Fibrocartilage + Inflammatory curated + top-15 data-driven markers each |
| Cartilage_context_specific | composite | Pre_hypertrophic curated + top-15 markers |
| Cystic_lesion_specific | composite | Hypertrophic_osteogenic curated + top-15 markers |

Data-driven markers: wilcox, min.pct=0.10, logfc.threshold=0.25, only.pos=TRUE, max.cells.per.ident=200, FDR<0.05, top 20 by avg_log2FC.

## 4. External validation cohorts

### GSE74089 (Agilent 4x44K, GPL13497; hip cartilage; NFH 4 vs Control 4)
- Probe→gene: Agilent Feature Extraction ProbeName→SystematicName → org.Hs.eg.db → symbol; 44,495 probes → 35,607 annotated → 20,401 unique genes; collapsing = highest-variance probe per gene
- DE: limma; signature testing: ssGSEA + mean z-score
- Re-run verification (revision): Spearman rho(logFC) = 1.0 vs original; 8,952/8,952 significant genes reproduced

### GSE345158 (Illumina NovaSeq 6000; femoral head bulk; SONFH 5 vs Control 4)
- Ensembl→symbol via org.Hs.eg.db; 18,834 genes; DESeq2; ssGSEA + mean z-score

## 5. Robustness analyses (Phase 1)

| Analysis | Parameters |
|----------|-----------|
| Leave-one-donor-out | 6 iterations, edgeR pseudobulk each |
| Cell-number sensitivity | Downsampling 25–90% of cells |
| Scoring sensitivity | AddModuleScore vs UCell concordance |
| Gene-dominance | Leave-one-gene-out per program (10k-cell subset) |
| Resolution sensitivity | Program scores across 7 resolutions (0.2–0.8) |

## 6. Revision-phase parameters (Phase 2)

### scDblFinder (R2 Task 3)
- Per-donor calling: `samples = donor_id`; GSE290411 SIONFH-4 excluded (n=1 cell < 10-cell minimum)
- Results: GSE311084 overall 9.94% (range 6.06–12.69%/donor); GSE290411 overall 4.40% (0–4.54%)

### Negative control (R2 Task 4)
- 1,000 random gene sets per program; matched to program size and expression level (deciles of log1p mean counts)
- UCell scoring on 10,000-cell stratified subset; two-sided empirical p = (rank + 1)/(N + 1)
- Comparison signatures: Hallmark EMT, generic ECM, generic fibroblast

### Meta-analysis (R1 #2)
- Model: random-effects, DerSimonian-Laird (metafor `rma(yi, sei, method="DL")`)
- Effect size: disease effect (ONFH − Control mean z-score difference) per external dataset
- 8 signatures × 2 datasets; heterogeneity: Q, I², τ²

### Donor-level statistics (R2 Task 1)
- Welch t-test (unequal variance) + Cohen's d, n=3 vs 3

### Cell composition (R2 Task 2)
- Per-sample raw/filtered counts; canonical marker pct-expressing within Fibrocartilage-dominant cells (SOX9, ACAN, COL2A1)

## 7. Reproducibility safeguards

1. Signatures frozen before validation data were examined (script + timestamped outputs)
2. Frozen signature file never altered post-freezing (215 genes, 8 signatures)
3. Donor-level inference throughout; cells never treated as independent replicates
4. Independent re-run of GSE74089 pipeline reproduced original exactly (rho = 1.0)
5. Negative-control null distributions archived (`R2_negative_control_null_distribution.csv`)
6. All figures exported as SVG (editable text) + PNG 300 dpi
7. `sessionInfo.txt` records exact R 4.4.2 package versions
