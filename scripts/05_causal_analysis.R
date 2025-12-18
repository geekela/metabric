
# 05_causal_analysis.R
# Analyze causal/directional edges from MIIC and compare across methods

library(tidyverse)

cat("=== Causal Edge Analysis ===\n\n")

# 1. Load All Results

miic_edges <- read_csv("results/miic_edges.csv", show_col_types = FALSE)
genie3_edges <- read_csv("results/genie3_edges.csv", show_col_types = FALSE)
aracne_edges <- read_csv("results/aracne_edges.csv", show_col_types = FALSE)

cat("Loaded edges:\n")
cat("  MIIC:", nrow(miic_edges), "\n")
cat("  GENIE3:", nrow(genie3_edges), "\n")
cat("  ARACNE:", nrow(aracne_edges), "\n")

# 2. Analyze MIIC Directional Evidence

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("MIIC CAUSAL/DIRECTIONAL ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Check available columns
cat("\nMIIC columns available:\n")
print(colnames(miic_edges))

# MIIC orientation columns:
# - p_x2y: probability of X -> Y direction
# - p_y2x: probability of Y -> X direction  
# - ort_inferred: orientation inferred (1 = yes)
# - is_causal: identified as causal

# Classify edges by direction strength
miic_directed <- miic_edges %>%
  mutate(
    # Direction strength: difference in directional probabilities
    direction_strength = abs(p_x2y - p_y2x),
    
    # Infer direction
    direction = case_when(
      p_x2y < 0.1 & p_y2x > 0.4 ~ paste(x, "->", y),
      p_y2x < 0.1 & p_x2y > 0.4 ~ paste(y, "->", x),
      p_x2y < 0.2 & p_y2x > 0.3 ~ paste(x, "->", y, "(weak)"),
      p_y2x < 0.2 & p_x2y > 0.3 ~ paste(y, "->", x, "(weak)"),
      TRUE ~ "undirected"
    ),
    
    # Source and target (for directed edges)
    source = case_when(
      p_x2y < p_y2x ~ x,
      p_y2x < p_x2y ~ y,
      TRUE ~ NA_character_
    ),
    target = case_when(
      p_x2y < p_y2x ~ y,
      p_y2x < p_x2y ~ x,
      TRUE ~ NA_character_
    ),
    
    # Is strongly directed?
    is_directed = (p_x2y < 0.2 | p_y2x < 0.2) & direction_strength > 0.2
  )

# Summary of directionality
cat("\nEdge directionality summary:\n")
direction_summary <- miic_directed %>%
  group_by(subtype) %>%
  summarise(
    total_edges = n(),
    directed_edges = sum(is_directed, na.rm = TRUE),
    pct_directed = round(100 * directed_edges / total_edges, 1),
    .groups = "drop"
  )
print(direction_summary)

# Strongly directed edges
strong_directed <- miic_directed %>%
  filter(is_directed) %>%
  select(subtype, edge_id, source, target, direction, p_x2y, p_y2x, direction_strength, info) %>%
  arrange(subtype, desc(direction_strength))

cat("\nStrongly directed edges per subtype:\n")
print(table(strong_directed$subtype))

# Save directed edges
write_csv(strong_directed, "results/miic_directed_edges.csv")

# 3. Identify Regulatory Hubs (Genes that regulate many others)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("REGULATORY HUB ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Count outgoing edges (gene as regulator)
regulator_counts <- strong_directed %>%
  filter(!is.na(source)) %>%
  group_by(subtype, source) %>%
  summarise(n_targets = n(), .groups = "drop") %>%
  arrange(subtype, desc(n_targets))

# Count incoming edges (gene as target)
target_counts <- strong_directed %>%
  filter(!is.na(target)) %>%
  group_by(subtype, target) %>%
  summarise(n_regulators = n(), .groups = "drop") %>%
  arrange(subtype, desc(n_regulators))

cat("\nTop regulators per subtype (genes with most outgoing edges):\n")
top_regulators <- regulator_counts %>%
  group_by(subtype) %>%
  slice_max(n_targets, n = 5) %>%
  ungroup()
print(top_regulators)

cat("\nTop targets per subtype (genes with most incoming edges):\n")
top_targets <- target_counts %>%
  group_by(subtype) %>%
  slice_max(n_regulators, n = 5) %>%
  ungroup()
print(top_targets)

# Save hub analysis
write_csv(regulator_counts, "results/miic_regulators.csv")
write_csv(target_counts, "results/miic_targets.csv")

# 4. Compare MIIC Directions with GENIE3

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("MIIC vs GENIE3 DIRECTION COMPARISON\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# GENIE3 is inherently directed: x (regulator) -> y (target)
# Compare with MIIC inferred directions

# Create standardized directed edge IDs
miic_dir_edges <- strong_directed %>%
  filter(!is.na(source), !is.na(target)) %>%
  mutate(directed_edge = paste(source, target, sep = "->"))

genie3_dir_edges <- genie3_edges %>%
  mutate(
    source = as.character(x),
    target = as.character(y),
    directed_edge = paste(source, target, sep = "->"),
    edge_id = paste(pmin(as.character(x), as.character(y)), 
                    pmax(as.character(x), as.character(y)), sep = "---")
  )

# Find agreements per subtype
direction_agreement <- data.frame()

for (st in unique(miic_dir_edges$subtype)) {
  miic_st <- miic_dir_edges %>% filter(subtype == st)
  genie3_st <- genie3_dir_edges %>% filter(subtype == st)
  
  # Same direction
  same_dir <- intersect(miic_st$directed_edge, genie3_st$directed_edge)
  
  # Opposite direction (MIIC says A->B, GENIE3 says B->A)
  miic_reversed <- paste(miic_st$target, miic_st$source, sep = "->")
  opposite_dir <- intersect(miic_reversed, genie3_st$directed_edge)
  
  # MIIC edges that overlap with GENIE3 (either direction)
  miic_undirected <- miic_st$edge_id
  genie3_undirected <- genie3_st$edge_id
  overlap <- intersect(miic_undirected, genie3_undirected)
  
  direction_agreement <- rbind(direction_agreement, data.frame(
    subtype = st,
    miic_directed = nrow(miic_st),
    same_direction = length(same_dir),
    opposite_direction = length(opposite_dir),
    edge_overlap = length(overlap)
  ))
}

cat("\nDirection agreement between MIIC and GENIE3:\n")
print(direction_agreement)

write_csv(direction_agreement, "results/direction_agreement.csv")

# 5. Three-Method Consensus (High-Confidence Edges)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("THREE-METHOD CONSENSUS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Find edges present in all 3 methods (per subtype)
consensus_results <- data.frame()

# Create proper edge_ids for GENIE3 and ARACNE
genie3_with_id <- genie3_edges %>%
  mutate(edge_id = paste(pmin(as.character(x), as.character(y)), 
                         pmax(as.character(x), as.character(y)), sep = "---"))

aracne_with_id <- aracne_edges %>%
  mutate(edge_id = paste(pmin(as.character(x), as.character(y)), 
                         pmax(as.character(x), as.character(y)), sep = "---"))

for (st in c("basal", "her2", "luma", "lumb")) {
  miic_st <- miic_edges %>% filter(subtype == st) %>% pull(edge_id) %>% unique()
  genie3_st <- genie3_with_id %>% filter(subtype == st) %>% pull(edge_id) %>% unique()
  aracne_st <- aracne_with_id %>% filter(subtype == st) %>% pull(edge_id) %>% unique()
  
  # Consensus
  two_method <- union(
    intersect(miic_st, genie3_st),
    union(intersect(miic_st, aracne_st), intersect(genie3_st, aracne_st))
  )
  three_method <- intersect(intersect(miic_st, genie3_st), aracne_st)
  
  consensus_results <- rbind(consensus_results, data.frame(
    subtype = st,
    miic = length(miic_st),
    genie3 = length(genie3_st),
    aracne = length(aracne_st),
    two_methods = length(two_method),
    three_methods = length(three_method)
  ))
  
  # Store 3-method consensus edges
  if (length(three_method) > 0) {
    consensus_edges <- miic_edges %>%
      filter(subtype == st, edge_id %in% three_method) %>%
      mutate(consensus = "3-method")
    
    if (st == "basal") {
      all_consensus <- consensus_edges
    } else {
      all_consensus <- rbind(all_consensus, consensus_edges)
    }
  }
}

cat("\nConsensus edge counts:\n")
print(consensus_results)

if (exists("all_consensus") && nrow(all_consensus) > 0) {
  write_csv(all_consensus, "results/consensus_3method_edges.csv")
  cat("\nSaved 3-method consensus edges:", nrow(all_consensus), "\n")
}

write_csv(consensus_results, "results/consensus_summary.csv")

# 6. Known Driver Gene Analysis (with direction)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("DRIVER GENE CAUSAL ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Known breast cancer drivers/regulators (lowercase to match data)
known_regulators <- tolower(c(
  # Transcription factors
  "ESR1", "GATA3", "FOXA1", "MYC", "E2F1", "E2F3", "E2F4", "RUNX1",
  # Signaling hubs  
  "TP53", "PIK3CA", "ERBB2", "AKT1", "PTEN", "RB1",
  # Kinases
  "CDK4", "CDK6", "MAPK1", "MAPK3", "JAK1", "JAK2", "SRC"
))

known_targets <- tolower(c(
  # Cell cycle
  "CCND1", "CCNE1", "CCNB1", "CDK2", "CDKN1A", "CDKN1B",
  # Apoptosis
  "BCL2L1", "CASP8", "CASP9",
  # Other downstream
  "MMP1", "MMP11", "VEGFA", "VEGFB"
))

# Check if known regulators are identified as sources
driver_as_regulator <- strong_directed %>%
  filter(tolower(source) %in% known_regulators) %>%
  select(subtype, source, target, direction_strength)

cat("\nKnown drivers acting as regulators (outgoing edges):\n")
if (nrow(driver_as_regulator) > 0) {
  print(driver_as_regulator %>% arrange(subtype, source))
} else {
  cat("None found in strongly directed edges\n")
}

# Check if known targets are identified as targets
known_as_target <- strong_directed %>%
  filter(tolower(target) %in% known_targets) %>%
  select(subtype, source, target, direction_strength)

cat("\nKnown targets receiving regulatory edges:\n")
if (nrow(known_as_target) > 0) {
  print(known_as_target %>% arrange(subtype, target))
} else {
  cat("None found in strongly directed edges\n")
}

# Broader check: any edges involving driver genes
all_drivers <- c(known_regulators, known_targets)
driver_edges <- miic_directed %>%
  filter(tolower(x) %in% all_drivers | tolower(y) %in% all_drivers) %>%
  select(subtype, x, y, direction, is_directed, info) %>%
  arrange(subtype, desc(info))

cat("\nAll edges involving known driver genes:\n")
cat("Total:", nrow(driver_edges), "\n")
print(head(driver_edges, 20))

write_csv(driver_edges, "results/driver_gene_edges.csv")

# 7. Subtype-Specific Causal Patterns

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUBTYPE-SPECIFIC CAUSAL PATTERNS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Find directed edges unique to each subtype
subtype_specific_directed <- strong_directed %>%
  group_by(edge_id) %>%
  mutate(n_subtypes = n_distinct(subtype)) %>%
  ungroup() %>%
  filter(n_subtypes == 1)

cat("\nSubtype-specific directed edges:\n")
print(table(subtype_specific_directed$subtype))

# Top subtype-specific regulatory relationships
cat("\nTop subtype-specific causal edges:\n")
for (st in c("basal", "her2", "luma", "lumb")) {
  cat("\n", toupper(st), ":\n")
  st_edges <- subtype_specific_directed %>%
    filter(subtype == st) %>%
    arrange(desc(direction_strength)) %>%
    head(5)
  if (nrow(st_edges) > 0) {
    print(st_edges %>% select(source, target, direction_strength))
  } else {
    cat("  No subtype-specific directed edges\n")
  }
}

write_csv(subtype_specific_directed, "results/subtype_specific_directed.csv")

# 8. Summary Report

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("\n1. MIIC DIRECTIONALITY:\n")
cat("   - Total MIIC edges:", nrow(miic_edges), "\n")
cat("   - Strongly directed:", nrow(strong_directed), 
    "(", round(100*nrow(strong_directed)/nrow(miic_edges), 1), "%)\n")

cat("\n2. METHOD AGREEMENT:\n")
cat("   - See consensus_summary.csv for details\n")

cat("\n3. KEY FINDINGS:\n")
cat("   - Top regulators and targets saved to results/\n")
cat("   - Driver gene edges saved to results/driver_gene_edges.csv\n")

cat("\n=== Causal Analysis Complete ===\n")
cat("\nFiles saved:\n")
cat("  - results/miic_directed_edges.csv\n")
cat("  - results/miic_regulators.csv\n")
cat("  - results/miic_targets.csv\n")
cat("  - results/direction_agreement.csv\n")
cat("  - results/consensus_summary.csv\n")
cat("  - results/consensus_3method_edges.csv\n")
cat("  - results/driver_gene_edges.csv\n")
cat("  - results/subtype_specific_directed.csv\n")
