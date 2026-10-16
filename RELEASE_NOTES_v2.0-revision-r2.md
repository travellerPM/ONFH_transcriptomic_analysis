# Release notes — v2.0-revision-r2 (Frontiers second revision, Reviewer 2)

**Source state at delivery: local Biomni version control only; this archive has not been uploaded to a remote repository.** A local git repository was initialized for the
analysis working copy (`/workspace/phase1`, commit `c7cf65b`, annotated tag
`v2.0-revision-r2`). **No push to the authors’ GitHub repository was performed as part of this delivery.** The commands below are prepared for
the authors to create the corresponding remote release.

## Contents of this release

1. Preserved original Phase-1 script: `scripts/07_freeze_signatures.R` (signature
   freezing before validation; frozen signatures unchanged).
2. Executed revision-phase scripts (Reviewer 1/2 analyses): `scripts/R1_*.R`,
   `scripts/R2_*.R`, incl. the corrected
   `R2_task3_doublet_cartilage_lineage_comparator.R` (dated correction note inside).
3. New integrity-audit scripts: `scripts/R2_task1_cell_annotation_provenance.R`,
   `scripts/R2_task1b_sensitivity_condition_composition.R`.
4. Partially RECONSTRUCTED Phase-1 modules (`scripts/reconstructed_phase1/`, all headed
   RECONSTRUCTED, none claimed to be the originals) with executed consistency
   checks (`RECONSTRUCTED_validation_check.R`; V1 AddModuleScore exact
   reproduction rho=1.000; V2 edgeR logFC rho=0.9999; V3-V5 structural checks).
5. Updated `README.md` (pipeline, provenance, annotation audit, validation table)
   and `parameters.md`.
6. Audit deliverables: `R2_cell_identity_audit.md`,
   `R2_final_integrity_audit.md`, `R2_Manuscript_Ready_Text.md`,
   `R2_Response_Evidence.md`, corrected statistics tables.

## Prepared commands (run by the authors with their own remote)

```bash
cd <analysis-working-copy>
# For existing repositories: clone your existing repository FIRST, copy upload files into it,
# commit + push that current branch; then create tag on the pushed commit.
# Never assume Biomni local commit c7cf65b exists in your existing GitHub repo.
# then create a GitHub release from v2.0-revision-r2. Attach selected sanitized
# results only if needed; do not automatically attach the entire Biomni archive.
```

## Environment

R 4.4.2; Seurat 5.5.1; GSVA (ssgseaParam, normalize=TRUE); edgeR; DESeq2; limma;
UCell; scDblFinder; metafor; seed 42 throughout. Full versions: `sessionInfo.txt`
and per-script sessionInfo files in `tables/`.

Note: Some RECONSTRUCTED modules contain workflow placeholders and unknown parameters. Do not claim complete end-to-end execution in the release description.
