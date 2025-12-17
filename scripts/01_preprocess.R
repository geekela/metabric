# =============================================================================
# 01_preprocess.R
# Load and preprocess METABRIC data for MIIC network inference
# =============================================================================

library(tidyverse)

cat("=== METABRIC Data Preprocessing ===\n\n")

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------

data_path <- "data/METABRIC_RNA_Mutation.csv"

if (!file.exists(data_path)) {
  stop("Data file not found. Please download from Kaggle and place in data/ folder.")
}

df <- read_csv(data_path, show_col_types = FALSE)
cat("Loaded data:", nrow(df), "samples x", ncol(df), "columns\n")

# -----------------------------------------------------------------------------
# 2. Identify Column Types
# -----------------------------------------------------------------------------

# Clinical columns (first ~31 columns based on Kaggle dataset structure)
clinical_cols <- c(
  "patient_id", "age_at_diagnosis", "type_of_breast_surgery", 
  "cancer_type", "cancer_type_detailed", "cellularity", "chemotherapy",
  "pam50_+_claudin-low_subtype", "cohort", "er_status_measured_by_ihc",
  "er_status", "neoplasm_histologic_grade", "her2_status_measured_by_snp6",
  "her2_status", "tumor_other_histologic_subtype", "hormone_therapy",
  "inferred_menopausal_state", "integrative_cluster", "primary_tumor_laterality",
  "lymph_nodes_examined_positive", "mutation_count", "nottingham_prognostic_index",
  "oncotree_code", "overall_survival_months", "overall_survival",
  "pr_status", "radio_therapy", "3-gene_classifier_subtype", "tumor_size",
  "tumor_stage", "death_from_cancer"
)

# Gene expression columns (exclude clinical and mutation columns)
all_cols <- colnames(df)
mutation_cols <- all_cols[grepl("_mut$", all_cols)]
gene_cols <- setdiff(all_cols, c(clinical_cols, mutation_cols))

cat("\nColumn breakdown:\n")
cat("  Clinical:", length(intersect(clinical_cols, all_cols)), "\n")
cat("  Gene expression:", length(gene_cols), "\n
")
cat("  Mutation:", length(mutation_cols), "\n")

# -----------------------------------------------------------------------------
# 3. Extract and Clean Expression Data
# -----------------------------------------------------------------------------

# Extract expression matrix
expr_df <- df %>% select(all_of(gene_cols))

# Check for non-numeric columns that slipped through
numeric_check <- sapply(expr_df, is.numeric)
if (!all(numeric_check)) {
  non_numeric <- names(numeric_check)[!numeric_check]
  cat("\nWarning: Removing non-numeric columns:", paste(non_numeric, collapse = ", "), "\n")
  expr_df <- expr_df %>% select(where(is.numeric))
  gene_cols <- colnames(expr_df)
}

cat("\nExpression matrix:", nrow(expr_df), "samples x", ncol(expr_df), "genes\n")

# -----------------------------------------------------------------------------
# 4. Handle Missing Values
# -----------------------------------------------------------------------------

# Check NA proportion per gene
na_prop <- colMeans(is.na(expr_df))
cat("\nNA proportion summary:\n")
print(summary(na_prop))

# Remove genes with >10% missing
genes_to_keep <- names(na_prop[na_prop < 0.1])
expr_df <- expr_df %>% select(all_of(genes_to_keep))
cat("\nGenes after NA filter:", ncol(expr_df), "\n")

# Impute remaining NAs with column median
expr_df <- expr_df %>%
  mutate(across(everything(), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

# -----------------------------------------------------------------------------
# 5. Variance Filtering
# -----------------------------------------------------------------------------

# Keep top N most variable genes (MIIC can't handle thousands efficiently)
N_GENES <- 300  # Adjust based on computational resources

gene_var <- apply(expr_df, 2, var, na.rm = TRUE)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(N_GENES, length(gene_var))]

expr_filtered <- expr_df %>% select(all_of(top_genes))
cat("\nAfter variance filter: top", ncol(expr_filtered), "genes\n")

# -----------------------------------------------------------------------------
# 6. Extract Clinical Info and Subtypes
# -----------------------------------------------------------------------------

# Get subtype column (handle different possible names)
subtype_col <- "pam50_+_claudin-low_subtype"
if (!subtype_col %in% colnames(df)) {
  subtype_col <- colnames(df)[grepl("pam50|claudin|subtype", colnames(df), ignore.case = TRUE)][1]
}

clinical_df <- df %>%
  select(patient_id, all_of(subtype_col), overall_survival_months, overall_survival) %>%
  rename(subtype = all_of(subtype_col))

cat("\nSubtype distribution:\n")
print(table(clinical_df$subtype, useNA = "ifany"))

# Map subtypes to cleaner names
subtype_mapping <- c(
  "LumA" = "LumA",
  "LumB" = "LumB", 
  "Her2" = "Her2",
  "Basal" = "Basal",
  "claudin-low" = "ClaudinLow",
  "Normal" = "Normal",
  "NC" = NA
)

clinical_df <- clinical_df %>%
  mutate(subtype_clean = case_when(
    grepl("LumA|Luminal A", subtype, ignore.case = TRUE) ~ "LumA",
    grepl("LumB|Luminal B", subtype, ignore.case = TRUE) ~ "LumB",
    grepl("Her2|HER2", subtype, ignore.case = TRUE) ~ "Her2",
    grepl("Basal", subtype, ignore.case = TRUE) ~ "Basal",
    grepl("claudin", subtype, ignore.case = TRUE) ~ "ClaudinLow",
    grepl("Normal", subtype, ignore.case = TRUE) ~ "Normal",
    TRUE ~ NA_character_
  ))

cat("\nCleaned subtype distribution:\n")
print(table(clinical_df$subtype_clean, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 7. Split Expression by Subtype
# -----------------------------------------------------------------------------

# Add patient_id to expression data
expr_with_id <- bind_cols(
  patient_id = df$patient_id,
  expr_filtered
)

# Merge with clinical
expr_with_subtype <- expr_with_id %>%
  left_join(clinical_df %>% select(patient_id, subtype_clean), by = "patient_id")

# Define subtypes to analyze (need sufficient sample size)
subtypes_to_analyze <- c("LumA", "LumB", "Her2", "Basal")

expr_by_subtype <- list()
for (st in subtypes_to_analyze) {
  expr_by_subtype[[st]] <- expr_with_subtype %>%
    filter(subtype_clean == st) %>%
    select(-patient_id, -subtype_clean)
  
  cat(st, ":", nrow(expr_by_subtype[[st]]), "samples\n")
}

# -----------------------------------------------------------------------------
# 8. Save Processed Data
# -----------------------------------------------------------------------------

# Save full filtered expression
saveRDS(expr_filtered, "data/expr_filtered.rds")

# Save expression by subtype
saveRDS(expr_by_subtype, "data/expr_by_subtype.rds")

# Save clinical data
saveRDS(clinical_df, "data/clinical_clean.rds")

# Save gene list
writeLines(top_genes, "data/selected_genes.txt")

cat("\n=== Preprocessing Complete ===\n")
cat("Saved:\n")
cat("  - data/expr_filtered.rds\n")
cat("  - data/expr_by_subtype.rds\n")
cat("  - data/clinical_clean.rds\n")
cat("  - data/selected_genes.txt\n")
