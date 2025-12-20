#!/usr/bin/env Rscript

# 00_feature_selection.R
# Improved gene selection with biological justification


library(tidyverse)

cat("=== Improved Feature Selection ===\n\n")


# 1. Load expression data

cat("Loading expression data...\n")

# Load the original METABRIC data
# Adjust path as needed
expr_data <- read.csv("data/METABRIC_RNA_Mutation.csv", stringsAsFactors = FALSE)

# Separate clinical and expression columns
clinical_cols <- c("patient_id", "age_at_diagnosis", "type_of_breast_surgery", 
                   "cancer_type", "cancer_type_detailed", "cellularity", 
                   "chemotherapy", "pam50_._claudin.low_subtype", "cohort",
                   "er_status_measured_by_ihc", "er_status", "neoplasm_histologic_grade",
                   "her2_status_measured_by_snp6", "her2_status", "tumor_other_histologic_subtype",
                   "hormone_therapy", "inferred_menopausal_state", "integrative_cluster",
                   "primary_tumor_laterality", "lymph_nodes_examined_positive",
                   "mutation_count", "nottingham_prognostic_index", "oncotree_code",
                   "overall_survival_months", "overall_survival", "pr_status",
                   "radio_therapy", "3.gene_classifier_subtype", "tumor_size", 
                   "tumor_stage", "death_from_cancer")

# Get expression columns (all others)
expr_cols <- setdiff(names(expr_data), clinical_cols)
cat(sprintf("Total genes in dataset: %d\n", length(expr_cols)))

# Extract expression matrix
expr_matrix <- expr_data[, expr_cols]
rownames(expr_matrix) <- expr_data$patient_id

# Get subtype information
subtype_col <- "pam50_._claudin.low_subtype"
subtypes <- expr_data[[subtype_col]]

# 2. Define biologically relevant gene sets

cat("\nDefining biological gene sets...\n")

# PAM50 signature genes (breast cancer subtype classifiers)
pam50_genes <- c(
  # Luminal genes
  "ESR1", "PGR", "BCL2", "GATA3", "FOXA1", "XBP1", "MLPH", "GPR160", "TMEM45B",
  # Proliferation genes
  "MKI67", "CCNB1", "CCNE1", "CDC20", "CDC6", "CDCA1", "CEP55", "EXO1", "KNTC2",
  "MELK", "NDC80", "NUF2", "ORC6L", "PTTG1", "RRM2", "TYMS", "UBE2C", "UBE2T",
  # HER2 genes
  "ERBB2", "GRB7", "FGFR4",
  # Basal genes
  "KRT5", "KRT14", "KRT17", "SFRP1", "CDH3", "EGFR",
  # Other markers
  "ACTR3B", "ANLN", "BAG1", "BIRC5", "BLVRA", "CENPF", "CXXC5", "MAPT", "MDM2",
  "MMP11", "MYBL2", "MYC", "NAT1", "PHGDH", "SLC39A6"
)

# Cancer Gene Census - top breast cancer drivers
# Source: COSMIC (https://cancer.sanger.ac.uk/census)
cancer_drivers <- c(
  # Major breast cancer drivers
  "TP53", "PIK3CA", "GATA3", "CDH1", "MAP3K1", "PTEN", "AKT1", "CBFB", "RUNX1",
  "TBX3", "CDKN1B", "RB1", "ERBB2", "MYC", "CCND1", "FGFR1", "MDM2", "BRCA1", "BRCA2",
  # Signaling pathway genes
  "KRAS", "HRAS", "NRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
  "EGFR", "ERBB3", "ERBB4", "IGF1R", "FGFR2", "FGFR3", "MET",
  # PI3K/AKT/mTOR pathway
  "PIK3R1", "AKT2", "AKT3", "MTOR", "TSC1", "TSC2", "RICTOR", "RPTOR",
  # Cell cycle
  "CDK4", "CDK6", "CDKN2A", "CDKN2B", "CCNE1", "CCND2", "CCND3", "E2F1", "E2F3",
  # Apoptosis
  "BCL2", "BCL2L1", "MCL1", "BAX", "BAK1", "CASP3", "CASP8", "CASP9",
  # DNA repair
  "ATM", "ATR", "CHEK1", "CHEK2", "RAD51", "PALB2",
  # Epigenetic regulators
  "KMT2C", "KMT2D", "ARID1A", "SETD2", "KDM6A", "DNMT3A",
  # JAK-STAT pathway
  "JAK1", "JAK2", "STAT3", "STAT5A", "STAT5B"
)

# Known transcription factors (subset relevant to breast cancer)
# These are the regulators we expect to see as SOURCE nodes
transcription_factors <- c(
  "ESR1", "PGR", "GATA3", "FOXA1", "MYC", "E2F1", "E2F3", "TP53", "STAT3",
  "STAT5A", "STAT5B", "RUNX1", "CBFB", "TBX3", "SOX4", "SOX9", "SOX10",
  "FOXO1", "FOXO3", "FOXM1", "ETS1", "ELF3", "GRHL2", "ZEB1", "ZEB2",
  "SNAI1", "SNAI2", "TWIST1", "TWIST2", "NOTCH1", "NOTCH2", "NOTCH3",
  "HIF1A", "ARNT", "NF1", "NFKB1", "RELA", "JUN", "FOS", "ATF3", "CREB1"
)

# Combine all biological gene sets
biological_genes <- unique(c(pam50_genes, cancer_drivers, transcription_factors))
cat(sprintf("Total biological gene set: %d genes\n", length(biological_genes)))

# 3. Apply selection criteria

cat("\nApplying selection criteria...\n")

# Convert column names to uppercase for matching
expr_cols_upper <- toupper(expr_cols)
biological_genes_upper <- toupper(biological_genes)

# Find overlap
biological_in_data <- expr_cols[expr_cols_upper %in% biological_genes_upper]
cat(sprintf("Biological genes found in data: %d / %d\n", 
            length(biological_in_data), length(biological_genes)))

# Calculate variance for all genes
gene_variance <- apply(expr_matrix, 2, var, na.rm = TRUE)

# Create selection summary
selection_df <- data.frame(
  gene = expr_cols,
  gene_upper = expr_cols_upper,
  variance = gene_variance,
  is_pam50 = expr_cols_upper %in% toupper(pam50_genes),
  is_driver = expr_cols_upper %in% toupper(cancer_drivers),
  is_tf = expr_cols_upper %in% toupper(transcription_factors),
  stringsAsFactors = FALSE
)

selection_df$is_biological <- selection_df$is_pam50 | selection_df$is_driver | selection_df$is_tf

# Selection strategy:
# 1. Keep ALL biological genes (regardless of variance)
# 2. Add top variance genes to reach target size

target_genes <- 300
n_biological <- sum(selection_df$is_biological)
n_additional <- max(0, target_genes - n_biological)

# Get biological genes
selected_biological <- selection_df$gene[selection_df$is_biological]

# Get top variance non-biological genes
non_biological <- selection_df %>%
  filter(!is_biological) %>%
  arrange(desc(variance)) %>%
  head(n_additional) %>%
  pull(gene)

# Final selection
selected_genes <- c(selected_biological, non_biological)

cat(sprintf("\nFinal selection: %d genes\n", length(selected_genes)))
cat(sprintf("  - Biological genes: %d\n", length(selected_biological)))
cat(sprintf("    - PAM50 markers: %d\n", sum(selection_df$is_pam50)))
cat(sprintf("    - Cancer drivers: %d\n", sum(selection_df$is_driver)))
cat(sprintf("    - Transcription factors: %d\n", sum(selection_df$is_tf)))
cat(sprintf("  - High-variance genes: %d\n", length(non_biological)))

# 4. Create selection justification table

cat("\nCreating selection justification...\n")

justification <- selection_df %>%
  filter(gene %in% selected_genes) %>%
  mutate(
    selection_reason = case_when(
      is_pam50 & is_driver & is_tf ~ "PAM50 + Driver + TF",
      is_pam50 & is_driver ~ "PAM50 + Driver",
      is_pam50 & is_tf ~ "PAM50 + TF",
      is_driver & is_tf ~ "Driver + TF",
      is_pam50 ~ "PAM50 signature",
      is_driver ~ "Cancer driver",
      is_tf ~ "Transcription factor",
      TRUE ~ "High variance"
    )
  ) %>%
  arrange(desc(is_biological), desc(variance))

# Summary by category
cat("\nSelection breakdown:\n")
print(table(justification$selection_reason))

# 5. Save outputs

cat("\nSaving outputs...\n")

# Save selected genes list
write.csv(justification, "data/selected_genes_justified.csv", row.names = FALSE)

# Save gene list only
writeLines(selected_genes, "data/selected_genes.txt")

# Create filtered expression matrices per subtype
subtypes_clean <- c("Basal", "Her2", "LumA", "LumB")
subtype_mapping <- c(
  "Basal" = "Basal",
  "claudin-low" = NA,  # Exclude
  "Her2" = "Her2", 
  "LumA" = "LumA",
  "LumB" = "LumB",
  "NC" = NA,  # Exclude
  "Normal" = NA  # Exclude
)

expr_by_subtype <- list()

for (st in subtypes_clean) {
  idx <- which(subtypes == st)
  if (length(idx) > 50) {  # Minimum sample size
    expr_by_subtype[[tolower(st)]] <- expr_matrix[idx, selected_genes]
    cat(sprintf("  %s: %d samples, %d genes\n", st, nrow(expr_by_subtype[[tolower(st)]]), ncol(expr_by_subtype[[tolower(st)]])))
  }
}

saveRDS(expr_by_subtype, "data/expr_by_subtype_justified.rds")

# 6. Print justification summary for report

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("FEATURE SELECTION JUSTIFICATION (for report)\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("
Our gene selection strategy combines biological relevance with statistical
criteria to ensure both interpretability and network quality:

1. PAM50 SIGNATURE GENES (n=", sum(selection_df$is_pam50), ")
   - Clinically validated markers that define breast cancer subtypes
   - Essential for subtype-specific network differences
   - Source: Parker et al., J Clin Oncol 2009

2. CANCER DRIVER GENES (n=", sum(selection_df$is_driver), ")
   - Known oncogenes and tumor suppressors from COSMIC Cancer Gene Census
   - Include major breast cancer drivers (TP53, PIK3CA, ERBB2, etc.)
   - Represent key signaling pathways (PI3K/AKT, MAPK, cell cycle)
   - Source: COSMIC v97

3. TRANSCRIPTION FACTORS (n=", sum(selection_df$is_tf), ")
   - Regulatory proteins that should appear as network hubs
   - Enable validation: TFs should have outgoing edges
   - Include breast cancer-specific TFs (ESR1, GATA3, FOXA1)

4. HIGH-VARIANCE GENES (n=", length(non_biological), ")
   - Additional genes showing high expression variability
   - Capture sample heterogeneity not covered by curated lists

TOTAL SELECTED: ", length(selected_genes), " genes

This approach ensures our network includes:
- Known regulatory relationships for validation
- Subtype-defining markers for biological interpretation
- Sufficient variability for statistical inference
", sep = "")

cat("\n=== Feature Selection Complete ===\n")
