#!/usr/bin/env R
# R1 #3: Tissue-context interpretation table
# Dataset × tissue × platform × context matrix for cross-cohort interpretation

cat("=== R1 #3: Tissue-Context Interpretation ===\n")

# --- Build dataset characteristics table ---
datasets <- data.frame(
  dataset = c("GSE311084", "GSE290411", "GSE74089", "GSE345158"),
  tissue = c("Femoral head cartilage", "Femoral head cartilage", "Hip articular cartilage", "Femoral head tissue"),
  tissue_context = c("Cartilage-enriched (chondrocytes sorted/selected)", "Cartilage-enriched (chondrocytes)", "Cartilage (articular surface)", "Bulk tissue (mixed: cartilage + bone + marrow)"),
  platform = c("Illumina NovaSeq 6000 (scRNA-seq)", "Illumina HiSeq 1000 (scRNA-seq)", "Agilent 4x44K (microarray)", "Illumina NovaSeq 6000 (bulk RNA-seq)"),
  data_type = c("Single-cell", "Single-cell", "Bulk microarray", "Bulk RNA-seq"),
  n_disease = c(3, 4, 4, 5),
  n_control = c(3, 0, 4, 4),
  disease_label = c("ONFH", "SIONFH", "NFH", "SONFH"),
  role = c("Discovery", "Internal cross-context", "External validation", "External validation"),
  stringsAsFactors = FALSE
)

cat("\n--- Dataset Characteristics ---\n")
print(datasets)

# --- Add lesion context ---
lesion_context <- data.frame(
  dataset = c("GSE311084", "GSE290411", "GSE74089", "GSE345158"),
  lesion_type = c("Necrotic + adjacent cartilage", "Necrotic femoral head", "Articular cartilage (NFH)", "Femoral head (SONFH)"),
  cystic_lesion = c("Not specified", "Not specified", "Not specified", "Not specified"),
  cartilage_specificity = c("High (chondrocyte compartment)", "High (chondrocyte compartment)", "High (cartilage tissue)", "Low (bulk tissue)"),
  expected_cell_composition = c("Chondrocytes + fibrocartilage candidates", "Chondrocytes (94.8% from 1 donor)", "Cartilage-resident cells", "Mixed: chondrocytes, osteoblasts, fibroblasts, immune, endothelial"),
  interpretation_caveat = c("Fibrocartilage program reflects chondrocyte-state remodeling", "Donor skew limits generalizability", "Microarray probes may not capture all genes", "Bulk signal diluted by non-chondrocyte cells"),
  stringsAsFactors = FALSE
)

cat("\n--- Lesion Context ---\n")
print(lesion_context)

# --- Add validation interpretation per dataset ---
val_interpretation <- data.frame(
  dataset = c("GSE74089", "GSE345158"),
  fibrocartilage_concordant = c("YES (p=0.0025)", "YES (p=0.25, direction concordant)"),
  fibrocartilage_per_gene_pct = c(95.7, 73.9),
  interpretation = c(
    "Strong support: hip cartilage shares chondrocyte lineage with femoral head; microarray detects ECM/collagen probes well",
    "Direction concordant but weaker magnitude: bulk femoral head tissue dilutes chondrocyte signal with bone/marrow/immune cells"
  ),
  tissue_context_explanation = c(
    "Cartilage-enriched tissue → chondrocyte signal preserved → strong concordance",
    "Bulk tissue → chondrocyte signal diluted → weaker but directionally consistent"
  ),
  stringsAsFactors = FALSE
)

cat("\n--- Validation Interpretation ---\n")
print(val_interpretation)

# --- Combine into one table ---
output <- list(
  dataset_characteristics = datasets,
  lesion_context = lesion_context,
  validation_interpretation = val_interpretation
)

# Save as CSV files
write.csv(datasets, "/mnt/results/phase1/tables/R1_tissue_context_datasets.csv", row.names = FALSE)
write.csv(lesion_context, "/mnt/results/phase1/tables/R1_tissue_context_lesion.csv", row.names = FALSE)
write.csv(val_interpretation, "/mnt/results/phase1/tables/R1_tissue_context_validation.csv", row.names = FALSE)

# Also save as a single Excel file
suppressMessages(library(openxlsx))
wb <- createWorkbook()
addWorksheet(wb, "Dataset_Characteristics")
writeData(wb, "Dataset_Characteristics", datasets)
addWorksheet(wb, "Lesion_Context")
writeData(wb, "Lesion_Context", lesion_context)
addWorksheet(wb, "Validation_Interpretation")
writeData(wb, "Validation_Interpretation", val_interpretation)
saveWorkbook(wb, "/workspace/R1_tissue_context_interpretation.xlsx", overwrite = TRUE)
system("cp /workspace/R1_tissue_context_interpretation.xlsx /mnt/results/phase1/R1_tissue_context_interpretation.xlsx")

cat("\nSaved R1_tissue_context_interpretation.xlsx\n")
cat("\n=== R1 #3 Complete ===\n")
