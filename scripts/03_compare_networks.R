
# 03_compare_networks.R
# Compare networks across methods (MIIC, GENIE3, ARACNE) and subtypes


library(tidyverse)
library(igraph)

cat("=== Multi-Method Network Comparison ===\n\n")


# 1. Load Results from All Methods


cat("Loading edge lists from all methods...\n")

# Load edges from each method
miic_edges <- tryCatch(
  read_csv("results/miic_edges.csv", show_col_types = FALSE),
  error = function(e) { cat("MIIC edges not found\n"); data.frame() }
)

genie3_edges <- tryCatch(
  read_csv("results/genie3_edges.csv", show_col_types = FALSE),
  error = function(e) { cat("GENIE3 edges not found\n"); data.frame() }
)

aracne_edges <- tryCatch(
  read_csv("results/aracne_edges.csv", show_col_types = FALSE),
  error = function(e) { cat("ARACNE edges not found\n"); data.frame() }
)

# Standardize column names and combine
standardize_edges <- function(edges, method_name) {
  if (nrow(edges) == 0) return(data.frame())
  
  # Ensure we have the required columns
  if (!"x" %in% colnames(edges) || !"y" %in% colnames(edges)) {
    return(data.frame())
  }
  
  edges %>%
    mutate(
      method = method_name,
      edge_id = paste(pmin(x, y), pmax(x, y), sep = "---")
    ) %>%
    select(x, y, edge_id, subtype, method, everything())
}

# Combine all edges
all_edges <- bind_rows(
  standardize_edges(miic_edges, "MIIC"),
  standardize_edges(genie3_edges, "GENIE3"),
  standardize_edges(aracne_edges, "ARACNE")
)

methods <- unique(all_edges$method)
subtypes <- unique(all_edges$subtype)

cat("\nMethods loaded:", paste(methods, collapse = ", "), "\n")
cat("Subtypes:", paste(subtypes, collapse = ", "), "\n")
cat("Total edges:", nrow(all_edges), "\n")

# Summary table
cat("\nEdges per method and subtype:\n")
print(table(all_edges$method, all_edges$subtype))

# Save combined edges
write_csv(all_edges, "results/all_edges_combined.csv")


# 2. Cross-Method Agreement (per subtype)


cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("CROSS-METHOD AGREEMENT\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# For each subtype, compute Jaccard similarity between methods
method_agreement <- data.frame()

for (st in subtypes) {
  cat("\nSubtype:", st, "\n")
  
  edges_by_method <- list()
  for (m in methods) {
    edges_by_method[[m]] <- all_edges %>%
      filter(subtype == st, method == m) %>%
      pull(edge_id) %>%
      unique()
  }
  
  # Pairwise Jaccard
  for (i in 1:(length(methods) - 1)) {
    for (j in (i + 1):length(methods)) {
      m1 <- methods[i]
      m2 <- methods[j]
      
      e1 <- edges_by_method[[m1]]
      e2 <- edges_by_method[[m2]]
      
      intersection <- length(intersect(e1, e2))
      union <- length(union(e1, e2))
      jaccard <- ifelse(union > 0, intersection / union, 0)
      
      method_agreement <- rbind(method_agreement, data.frame(
        subtype = st,
        method1 = m1,
        method2 = m2,
        edges_m1 = length(e1),
        edges_m2 = length(e2),
        intersection = intersection,
        union = union,
        jaccard = jaccard
      ))
      
      cat(sprintf("  %s vs %s: Jaccard = %.3f (overlap: %d edges)\n", 
                  m1, m2, jaccard, intersection))
    }
  }
}

write_csv(method_agreement, "results/method_agreement.csv")


# 3. Consensus Edges (found by multiple methods)


cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("CONSENSUS EDGES\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Count methods per edge (within each subtype)
edge_method_count <- all_edges %>%
  group_by(subtype, edge_id) %>%
  summarise(
    n_methods = n_distinct(method),
    methods = paste(sort(unique(method)), collapse = ","),
    gene1 = first(x),
    gene2 = first(y),
    .groups = "drop"
  )

# Consensus edges (found by 2+ methods)
consensus_edges <- edge_method_count %>%
  filter(n_methods >= 2)

cat("\nConsensus edges (found by 2+ methods):\n")
print(table(consensus_edges$subtype, consensus_edges$n_methods))

# High-confidence edges (found by all 3 methods)
high_conf_edges <- edge_method_count %>%
  filter(n_methods == length(methods))

cat("\nHigh-confidence edges (all", length(methods), "methods):\n")
print(table(high_conf_edges$subtype))

# Save
write_csv(edge_method_count, "results/edge_method_count.csv")
write_csv(consensus_edges, "results/consensus_edges.csv")


# 4. Cross-Subtype Conservation (per method)


cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("CROSS-SUBTYPE CONSERVATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

subtype_conservation <- data.frame()

for (m in methods) {
  cat("\nMethod:", m, "\n")
  
  # Count subtypes per edge
  edge_subtype_count <- all_edges %>%
    filter(method == m) %>%
    group_by(edge_id) %>%
    summarise(
      n_subtypes = n_distinct(subtype),
      subtypes = paste(sort(unique(subtype)), collapse = ","),
      gene1 = first(x),
      gene2 = first(y),
      .groups = "drop"
    )
  
  # Conservation summary
  conservation_summary <- edge_subtype_count %>%
    count(n_subtypes) %>%
    mutate(
      method = m,
      category = case_when(
        n_subtypes == length(subtypes) ~ "Universal",
        n_subtypes >= 3 ~ "Conserved (3+)",
        n_subtypes == 2 ~ "Partial (2)",
        n_subtypes == 1 ~ "Subtype-specific"
      )
    )
  
  print(conservation_summary)
  subtype_conservation <- rbind(subtype_conservation, conservation_summary)
}

write_csv(subtype_conservation, "results/subtype_conservation.csv")


# 5. Hub Gene Analysis (per method and subtype)


cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("HUB GENE ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

compute_degree <- function(edges_df) {
  if (nrow(edges_df) == 0) return(data.frame(gene = character(), degree = integer()))
  
  g <- graph_from_data_frame(edges_df %>% select(x, y), directed = FALSE)
  g <- simplify(g)
  
  data.frame(
    gene = V(g)$name,
    degree = degree(g)
  )
}

hub_metrics <- data.frame()

for (m in methods) {
  for (st in subtypes) {
    edges_ms <- all_edges %>% filter(method == m, subtype == st)
    degrees <- compute_degree(edges_ms)
    
    if (nrow(degrees) > 0) {
      degrees$method <- m
      degrees$subtype <- st
      hub_metrics <- rbind(hub_metrics, degrees)
    }
  }
}

# Pivot for comparison
hub_wide <- hub_metrics %>%
  mutate(method_subtype = paste(method, subtype, sep = "_")) %>%
  select(gene, degree, method_subtype) %>%
  pivot_wider(names_from = method_subtype, values_from = degree, values_fill = 0)

# Compute stats
hub_wide <- hub_wide %>%
  rowwise() %>%
  mutate(
    total_degree = sum(c_across(-gene)),
    max_degree = max(c_across(-gene)),
    n_networks = sum(c_across(-gene) > 0)
  ) %>%
  ungroup() %>%
  arrange(desc(total_degree))

cat("\nTop 20 hub genes across all methods/subtypes:\n")
print(hub_wide %>% head(20) %>% select(gene, total_degree, max_degree, n_networks))

# Save
write_csv(hub_metrics, "results/hub_metrics_all.csv")
write_csv(hub_wide, "results/hub_comparison.csv")


# 6. Method Agreement on Hub Genes

cat("\n--- Hub Agreement Across Methods ---\n")

# Top 20 hubs per method (across all subtypes)
top_hubs_by_method <- hub_metrics %>%
  group_by(method, gene) %>%
  summarise(total_degree = sum(degree), .groups = "drop") %>%
  group_by(method) %>%
  slice_max(total_degree, n = 20) %>%
  ungroup()

# Find common hubs
hub_overlap <- top_hubs_by_method %>%
  group_by(gene) %>%
  summarise(
    n_methods = n_distinct(method),
    methods = paste(sort(unique(method)), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(n_methods >= 2) %>%
  arrange(desc(n_methods))

cat("\nGenes in top-20 hubs of multiple methods:\n")
print(hub_overlap)

write_csv(hub_overlap, "results/hub_overlap.csv")

# 7. Known Driver Gene Analysis

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("KNOWN DRIVER GENE ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

known_drivers <- c(
  "TP53", "PIK3CA", "ERBB2", "ESR1", "MYC", "CCND1", "GATA3", "CDH1",
  "PTEN", "RB1", "BRCA1", "BRCA2", "MAP3K1", "CDK4", "CDK6", "FGFR1",
  "AKT1", "FOXA1", "CBFB", "TBX3", "RUNX1", "NF1", "AURKA", "E2F1"
)

# Check which drivers are in our networks
genes_in_network <- unique(c(all_edges$x, all_edges$y))
drivers_found <- intersect(known_drivers, genes_in_network)

cat("Known drivers found in networks:", length(drivers_found), "/", length(known_drivers), "\n")
cat("Found:", paste(drivers_found, collapse = ", "), "\n")

# Driver degrees by method and subtype
if (length(drivers_found) > 0) {
  driver_analysis <- hub_metrics %>%
    filter(gene %in% drivers_found)
  
  driver_summary <- driver_analysis %>%
    pivot_wider(
      id_cols = gene,
      names_from = c(method, subtype),
      values_from = degree,
      values_fill = 0
    )
  
  cat("\nDriver gene degrees:\n")
  print(driver_summary)
  
  write_csv(driver_summary, "results/driver_gene_analysis.csv")
}

# 8. Network Statistics Summary

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("NETWORK STATISTICS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

network_stats <- data.frame()

for (m in methods) {
  for (st in subtypes) {
    edges_ms <- all_edges %>% filter(method == m, subtype == st)
    
    if (nrow(edges_ms) > 0) {
      g <- graph_from_data_frame(edges_ms %>% select(x, y), directed = FALSE)
      g <- simplify(g)
      
      network_stats <- rbind(network_stats, data.frame(
        method = m,
        subtype = st,
        n_nodes = vcount(g),
        n_edges = ecount(g),
        density = edge_density(g),
        avg_degree = mean(degree(g)),
        clustering = transitivity(g, type = "global"),
        n_components = components(g)$no
      ))
    }
  }
}

cat("\nNetwork statistics:\n")
print(network_stats)

write_csv(network_stats, "results/network_stats.csv")


# 9. Summary Report


cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY REPORT\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("\n1. METHODS COMPARED:", paste(methods, collapse = ", "), "\n")

cat("\n2. EDGES PER METHOD:\n")
for (m in methods) {
  n <- sum(all_edges$method == m)
  cat(sprintf("   %s: %d edges\n", m, n))
}

cat("\n3. CROSS-METHOD AGREEMENT (avg Jaccard per subtype):\n")
avg_jaccard <- method_agreement %>%
  group_by(subtype) %>%
  summarise(avg_jaccard = mean(jaccard))
print(avg_jaccard)

cat("\n4. CONSENSUS EDGES:\n")
cat(sprintf("   2+ methods: %d edges\n", sum(edge_method_count$n_methods >= 2)))
cat(sprintf("   All methods: %d edges\n", sum(edge_method_count$n_methods == length(methods))))

cat("\n5. TOP DIFFERENTIAL HUBS:\n")
top_5 <- hub_wide %>% head(5)
for (i in 1:nrow(top_5)) {
  cat(sprintf("   %s (total degree: %d)\n", top_5$gene[i], top_5$total_degree[i]))
}

cat("\n=== Comparison Analysis Complete ===\n")
cat("Results saved in results/ folder\n")
