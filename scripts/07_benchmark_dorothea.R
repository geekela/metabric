#!/usr/bin/env Rscript

# 07_benchmark_dorothea.R
# Quantitative benchmarking against DoRothEA reference network


library(tidyverse)

cat("=== Benchmarking Against DoRothEA ===\n\n")


# Helper function to standardize column names

standardize_edge_cols <- function(df) {
  cols <- tolower(names(df))
  names(df) <- cols
  
  # Rename source column
  if ("x" %in% cols && !"source" %in% cols) {
    names(df)[names(df) == "x"] <- "source"
  } else if ("from" %in% cols && !"source" %in% cols) {
    names(df)[names(df) == "from"] <- "source"
  } else if ("tf" %in% cols && !"source" %in% cols) {
    names(df)[names(df) == "tf"] <- "source"
  } else if ("regulator" %in% cols && !"source" %in% cols) {
    names(df)[names(df) == "regulator"] <- "source"
  }
  
  # Rename target column
  if ("y" %in% cols && !"target" %in% cols) {
    names(df)[names(df) == "y"] <- "target"
  } else if ("to" %in% cols && !"target" %in% cols) {
    names(df)[names(df) == "to"] <- "target"
  }
  
  return(df)
}


# 1. Install/Load DoRothEA

cat("Loading DoRothEA reference database...\n")

if (!requireNamespace("dorothea", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("dorothea")
}

library(dorothea)

data("dorothea_hs", package = "dorothea")

reference_edges <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  mutate(
    source_ref = toupper(tf),
    target_ref = toupper(target),
    edge_id = paste(pmin(source_ref, target_ref), pmax(source_ref, target_ref), sep = "_"),
    directed_id = paste(source_ref, target_ref, sep = "->")
  )

cat(sprintf("DoRothEA reference edges (confidence A-C): %d\n", nrow(reference_edges)))
cat(sprintf("  - Confidence A: %d\n", sum(reference_edges$confidence == "A")))
cat(sprintf("  - Confidence B: %d\n", sum(reference_edges$confidence == "B")))
cat(sprintf("  - Confidence C: %d\n", sum(reference_edges$confidence == "C")))
cat(sprintf("  - Unique TFs: %d\n", n_distinct(reference_edges$source_ref)))
cat(sprintf("  - Unique targets: %d\n", n_distinct(reference_edges$target_ref)))

# 2. Load inferred networks

cat("\nLoading inferred networks...\n")

# MIIC
miic_edges <- read.csv("results/miic_edges.csv", stringsAsFactors = FALSE)
cat(sprintf("  MIIC raw columns: %s\n", paste(names(miic_edges), collapse = ", ")))
miic_edges <- standardize_edge_cols(miic_edges)
miic_edges <- miic_edges %>%
  mutate(
    source = toupper(source),
    target = toupper(target),
    edge_id = paste(pmin(source, target), pmax(source, target), sep = "_")
  )

# GENIE3
genie3_edges <- read.csv("results/genie3_edges.csv", stringsAsFactors = FALSE)
cat(sprintf("  GENIE3 raw columns: %s\n", paste(names(genie3_edges), collapse = ", ")))
genie3_edges <- standardize_edge_cols(genie3_edges)
genie3_edges <- genie3_edges %>%
  mutate(
    source = toupper(source),
    target = toupper(target),
    edge_id = paste(pmin(source, target), pmax(source, target), sep = "_"),
    directed_id = paste(source, target, sep = "->")
  )

# ARACNE
aracne_edges <- read.csv("results/aracne_edges.csv", stringsAsFactors = FALSE)
cat(sprintf("  ARACNE raw columns: %s\n", paste(names(aracne_edges), collapse = ", ")))
aracne_edges <- standardize_edge_cols(aracne_edges)
aracne_edges <- aracne_edges %>%
  mutate(
    source = toupper(source),
    target = toupper(target),
    edge_id = paste(pmin(source, target), pmax(source, target), sep = "_")
  )

# Optional: directed MIIC edges
miic_directed <- NULL
if (file.exists("results/miic_directed_edges.csv")) {
  miic_directed <- read.csv("results/miic_directed_edges.csv", stringsAsFactors = FALSE)
  miic_directed <- standardize_edge_cols(miic_directed)
  miic_directed <- miic_directed %>%
    mutate(
      source = toupper(source),
      target = toupper(target),
      directed_id = paste(source, target, sep = "->")
    )
}

# Optional: consensus edges
consensus_edges <- NULL
if (file.exists("results/consensus_3method_edges.csv")) {
  consensus_edges <- read.csv("results/consensus_3method_edges.csv", stringsAsFactors = FALSE)
  consensus_edges <- standardize_edge_cols(consensus_edges)
  consensus_edges <- consensus_edges %>%
    mutate(
      source = toupper(source),
      target = toupper(target),
      edge_id = paste(pmin(source, target), pmax(source, target), sep = "_")
    )
}

cat(sprintf("\nLoaded edges:\n"))
cat(sprintf("  MIIC: %d edges\n", nrow(miic_edges)))
cat(sprintf("  GENIE3: %d edges\n", nrow(genie3_edges)))
cat(sprintf("  ARACNE: %d edges\n", nrow(aracne_edges)))
if (!is.null(miic_directed)) cat(sprintf("  MIIC directed: %d edges\n", nrow(miic_directed)))
if (!is.null(consensus_edges)) cat(sprintf("  3-method consensus: %d edges\n", nrow(consensus_edges)))

# 3. Benchmarking functions

benchmark_undirected <- function(inferred_edges, reference_edges, name) {
  inferred_ids <- unique(inferred_edges$edge_id)
  reference_ids <- unique(reference_edges$edge_id)
  
  tp <- length(intersect(inferred_ids, reference_ids))
  fp <- length(setdiff(inferred_ids, reference_ids))
  
  inferred_genes <- unique(c(inferred_edges$source, inferred_edges$target))
  reference_in_scope <- reference_edges %>%
    filter(source_ref %in% inferred_genes & target_ref %in% inferred_genes)
  reference_ids_in_scope <- unique(reference_in_scope$edge_id)
  fn <- length(setdiff(reference_ids_in_scope, inferred_ids))
  
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  f1 <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  
  data.frame(
    method = name, comparison = "undirected",
    inferred_edges = length(inferred_ids), reference_in_scope = length(reference_ids_in_scope),
    TP = tp, FP = fp, FN = fn,
    precision = round(precision, 4), recall = round(recall, 4), F1 = round(f1, 4)
  )
}

benchmark_directed <- function(inferred_edges, reference_edges, name) {
  inferred_ids <- unique(inferred_edges$directed_id)
  reference_ids <- unique(reference_edges$directed_id)
  
  tp <- length(intersect(inferred_ids, reference_ids))
  fp <- length(setdiff(inferred_ids, reference_ids))
  
  inferred_genes <- unique(c(inferred_edges$source, inferred_edges$target))
  reference_in_scope <- reference_edges %>%
    filter(source_ref %in% inferred_genes & target_ref %in% inferred_genes)
  reference_ids_in_scope <- unique(reference_in_scope$directed_id)
  fn <- length(setdiff(reference_ids_in_scope, inferred_ids))
  
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  f1 <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  
  data.frame(
    method = name, comparison = "directed",
    inferred_edges = length(inferred_ids), reference_in_scope = length(reference_ids_in_scope),
    TP = tp, FP = fp, FN = fn,
    precision = round(precision, 4), recall = round(recall, 4), F1 = round(f1, 4)
  )
}

# 4. Run benchmarking

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("BENCHMARKING RESULTS\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

results <- list()
results$miic <- benchmark_undirected(miic_edges, reference_edges, "MIIC")
results$genie3 <- benchmark_undirected(genie3_edges, reference_edges, "GENIE3")
results$aracne <- benchmark_undirected(aracne_edges, reference_edges, "ARACNE")
results$genie3_dir <- benchmark_directed(genie3_edges, reference_edges, "GENIE3 (directed)")

if (!is.null(miic_directed) && nrow(miic_directed) > 0) {
  results$miic_dir <- benchmark_directed(miic_directed, reference_edges, "MIIC (causal)")
}
if (!is.null(consensus_edges) && nrow(consensus_edges) > 0) {
  results$consensus <- benchmark_undirected(consensus_edges, reference_edges, "3-Method Consensus")
}

benchmark_df <- bind_rows(results)

cat("UNDIRECTED (edge exists regardless of direction):\n\n")
print(benchmark_df %>% filter(comparison == "undirected") %>% 
        select(method, inferred_edges, TP, precision, recall, F1) %>%
        as.data.frame(), row.names = FALSE)

cat("\nDIRECTED (TF->target must match):\n\n")
dir_res <- benchmark_df %>% filter(comparison == "directed")
if (nrow(dir_res) > 0) {
  print(dir_res %>% select(method, inferred_edges, TP, precision, recall, F1) %>%
          as.data.frame(), row.names = FALSE)
}

# 5. Per-subtype analysis

if ("subtype" %in% names(miic_edges)) {
  cat("\n", paste(rep("-", 50), collapse = ""), "\n")
  cat("PER-SUBTYPE RESULTS\n")
  cat(paste(rep("-", 50), collapse = ""), "\n\n")
  
  subtypes <- unique(miic_edges$subtype)
  subtype_results <- list()
  
  for (st in subtypes) {
    subtype_results[[paste0("miic_", st)]] <- benchmark_undirected(
      miic_edges %>% filter(subtype == st), reference_edges, paste("MIIC", st))
    subtype_results[[paste0("genie3_", st)]] <- benchmark_undirected(
      genie3_edges %>% filter(subtype == st), reference_edges, paste("GENIE3", st))
    subtype_results[[paste0("aracne_", st)]] <- benchmark_undirected(
      aracne_edges %>% filter(subtype == st), reference_edges, paste("ARACNE", st))
  }
  
  subtype_df <- bind_rows(subtype_results)
  print(subtype_df %>% select(method, TP, precision, recall) %>% as.data.frame(), row.names = FALSE)
  write.csv(subtype_df, "results/benchmark_per_subtype.csv", row.names = FALSE)
}

# 6. Show validated edges

cat("\n", paste(rep("-", 50), collapse = ""), "\n")
cat("VALIDATED EDGES (in DoRothEA)\n")
cat(paste(rep("-", 50), collapse = ""), "\n\n")

validated_miic <- miic_edges %>%
  filter(edge_id %in% reference_edges$edge_id) %>%
  left_join(reference_edges %>% select(edge_id, confidence) %>% distinct(), by = "edge_id")

cat(sprintf("MIIC: %d validated (%.1f%%)\n", nrow(validated_miic), 
            100 * nrow(validated_miic) / nrow(miic_edges)))

validated_genie3 <- genie3_edges %>% filter(edge_id %in% reference_edges$edge_id)
cat(sprintf("GENIE3: %d validated (%.1f%%)\n", nrow(validated_genie3),
            100 * nrow(validated_genie3) / nrow(genie3_edges)))

validated_aracne <- aracne_edges %>% filter(edge_id %in% reference_edges$edge_id)
cat(sprintf("ARACNE: %d validated (%.1f%%)\n", nrow(validated_aracne),
            100 * nrow(validated_aracne) / nrow(aracne_edges)))

if (!is.null(consensus_edges)) {
  validated_consensus <- consensus_edges %>% filter(edge_id %in% reference_edges$edge_id)
  cat(sprintf("Consensus: %d validated (%.1f%%)\n", nrow(validated_consensus),
              100 * nrow(validated_consensus) / nrow(consensus_edges)))
}

# Show examples
if (nrow(validated_miic) > 0) {
  cat("\nExample validated edges (MIIC, confidence A):\n")
  examples <- validated_miic %>% 
    filter(confidence == "A") %>%
    head(10) %>%
    select(source, target, confidence, any_of(c("subtype", "info")))
  print(as.data.frame(examples), row.names = FALSE)
}

# 7. Save results

cat("\nSaving results...\n")
write.csv(benchmark_df, "results/benchmark_overall.csv", row.names = FALSE)
write.csv(validated_miic, "results/validated_miic_edges.csv", row.names = FALSE)

# 8. Create figure

if (!dir.exists("figures")) dir.create("figures")

plot_data <- benchmark_df %>%
  filter(comparison == "undirected") %>%
  select(method, precision, recall, F1) %>%
  pivot_longer(cols = c(precision, recall, F1), names_to = "metric", values_to = "value")

p <- ggplot(plot_data, aes(x = reorder(method, -value), y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value)), 
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("precision" = "#2E86AB", "recall" = "#A23B72", "F1" = "#F18F01")) +
  labs(title = "Benchmark vs DoRothEA (confidence A-C)",
       subtitle = "Undirected comparison", x = "", y = "Score", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylim(0, max(plot_data$value) * 1.25)

ggsave("figures/benchmark_dorothea.pdf", p, width = 9, height = 5)
ggsave("figures/benchmark_dorothea.png", p, width = 9, height = 5, dpi = 150)
cat("Figure saved: figures/benchmark_dorothea.pdf\n")

# Summary

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY FOR REPORT\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(sprintf("
Reference: DoRothEA v1.20 (Garcia-Alonso et al., Genome Res 2019)
  - %d TF-target edges (confidence A-C)
  - %d unique transcription factors

Best precision: %s (%.1f%%)
Most validated edges: %s (%d TP)

Interpretation:
- Precision = %% of our edges found in DoRothEA
- Low recall expected (we use ~300 genes, DoRothEA has ~5000 targets)
- Higher precision for consensus = multi-method agreement is more reliable
",
nrow(reference_edges), n_distinct(reference_edges$source_ref),
benchmark_df$method[which.max(benchmark_df$precision)],
max(benchmark_df$precision) * 100,
benchmark_df$method[which.max(benchmark_df$TP)],
max(benchmark_df$TP)))

cat("\n=== Done ===\n")
