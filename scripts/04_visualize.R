
# 04_visualize.R
# Generate visualizations for multi-method network comparison


library(tidyverse)
library(igraph)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

cat("=== Generating Visualizations ===\n\n")

# Set theme
theme_set(theme_minimal(base_size = 12))

# Color palettes
method_colors <- c("MIIC" = "#E41A1C", "GENIE3" = "#377EB8", "ARACNE" = "#4DAF4A")
subtype_colors <- c("LumA" = "#1B9E77", "LumB" = "#D95F02", 
                    "Her2" = "#7570B3", "Basal" = "#E7298A")

# 1. Load Results

all_edges <- read_csv("results/all_edges_combined.csv", show_col_types = FALSE)
method_agreement <- read_csv("results/method_agreement.csv", show_col_types = FALSE)
consensus_edges <- read_csv("results/consensus_edges.csv", show_col_types = FALSE)
hub_metrics <- read_csv("results/hub_metrics_all.csv", show_col_types = FALSE)
hub_comparison <- read_csv("results/hub_comparison.csv", show_col_types = FALSE)
network_stats <- read_csv("results/network_stats.csv", show_col_types = FALSE)

methods <- unique(all_edges$method)
subtypes <- unique(all_edges$subtype)

# 2. Network Size Comparison (Method x Subtype)

cat("Creating network size comparison...\n")

p1 <- network_stats %>%
  ggplot(aes(x = subtype, y = n_edges, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Network Size by Method and Subtype",
    subtitle = "Number of edges inferred",
    x = "Subtype",
    y = "Number of Edges",
    fill = "Method"
  ) +
  theme(legend.position = "bottom")

ggsave("figures/01_network_size_comparison.pdf", p1, width = 10, height = 6)
ggsave("figures/01_network_size_comparison.png", p1, width = 10, height = 6, dpi = 150)

# 3. Method Agreement Heatmap

cat("Creating method agreement heatmap...\n")

# Create matrix for heatmap
agreement_matrix <- method_agreement %>%
  mutate(pair = paste(method1, method2, sep = " vs ")) %>%
  select(subtype, pair, jaccard) %>%
  pivot_wider(names_from = pair, values_from = jaccard) %>%
  column_to_rownames("subtype") %>%
  as.matrix()

pdf("figures/02_method_agreement_heatmap.pdf", width = 8, height = 5)
pheatmap(
  agreement_matrix,
  main = "Cross-Method Agreement (Jaccard Index)",
  color = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
  display_numbers = TRUE,
  number_format = "%.3f",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 12
)
dev.off()

png("figures/02_method_agreement_heatmap.png", width = 800, height = 500)
pheatmap(
  agreement_matrix,
  main = "Cross-Method Agreement (Jaccard Index)",
  color = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
  display_numbers = TRUE,
  number_format = "%.3f",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 12
)
dev.off()

# 4. Consensus Edges Bar Plot

cat("Creating consensus edges plot...\n")

consensus_summary <- consensus_edges %>%
  count(subtype, n_methods) %>%
  mutate(n_methods = factor(n_methods, labels = c("2 methods", "3 methods")[1:n_distinct(n_methods)]))

p2 <- ggplot(consensus_summary, aes(x = subtype, y = n, fill = n_methods)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Blues") +
  labs(
    title = "Consensus Edges Across Methods",
    subtitle = "Edges found by multiple inference methods",
    x = "Subtype",
    y = "Number of Edges",
    fill = "Agreement"
  ) +
  theme(legend.position = "bottom")

ggsave("figures/03_consensus_edges.pdf", p2, width = 8, height = 5)
ggsave("figures/03_consensus_edges.png", p2, width = 8, height = 5, dpi = 150)

# 5. Hub Gene Comparison Heatmap

cat("Creating hub gene heatmap...\n")

# Top 30 hub genes
top_hubs <- hub_comparison %>%
  head(30)

# Select only method_subtype columns for heatmap
hub_cols <- colnames(top_hubs)[grepl("_", colnames(top_hubs)) & 
                                 !colnames(top_hubs) %in% c("total_degree", "max_degree", "n_networks")]

hub_matrix <- top_hubs %>%
  select(gene, all_of(hub_cols)) %>%
  column_to_rownames("gene") %>%
  as.matrix()

pdf("figures/04_hub_gene_heatmap.pdf", width = 14, height = 10)
pheatmap(
  hub_matrix,
  main = "Top 30 Hub Genes: Degree Across Methods and Subtypes",
  color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 8
)
dev.off()

png("figures/04_hub_gene_heatmap.png", width = 1400, height = 1000)
pheatmap(
  hub_matrix,
  main = "Top 30 Hub Genes: Degree Across Methods and Subtypes",
  color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 8
)
dev.off()

# 6. Degree Distribution by Method

cat("Creating degree distribution plots...\n")

p3 <- hub_metrics %>%
  ggplot(aes(x = degree, fill = method)) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.6) +
  facet_grid(subtype ~ method) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Degree Distribution by Method and Subtype",
    x = "Node Degree",
    y = "Frequency"
  ) +
  theme(legend.position = "none")

ggsave("figures/05_degree_distribution.pdf", p3, width = 12, height = 10)
ggsave("figures/05_degree_distribution.png", p3, width = 12, height = 10, dpi = 150)

# 7. Average Degree Comparison

cat("Creating average degree comparison...\n")

p4 <- network_stats %>%
  ggplot(aes(x = method, y = avg_degree, fill = subtype)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = subtype_colors) +
  labs(
    title = "Average Node Degree by Method and Subtype",
    x = "Method",
    y = "Average Degree",
    fill = "Subtype"
  ) +
  theme(legend.position = "bottom")

ggsave("figures/06_avg_degree_comparison.pdf", p4, width = 10, height = 6)
ggsave("figures/06_avg_degree_comparison.png", p4, width = 10, height = 6, dpi = 150)

# 8. Method Comparison Radar/Summary

cat("Creating method summary plot...\n")

method_summary <- network_stats %>%
  group_by(method) %>%
  summarise(
    total_edges = sum(n_edges),
    avg_nodes = mean(n_nodes),
    avg_density = mean(density),
    avg_clustering = mean(clustering, na.rm = TRUE),
    .groups = "drop"
  )

p5 <- method_summary %>%
  pivot_longer(-method, names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = metric, y = value, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = method_colors) +
  facet_wrap(~metric, scales = "free") +
  labs(
    title = "Method Comparison: Network Properties",
    x = "",
    y = "Value",
    fill = "Method"
  ) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom")

ggsave("figures/07_method_summary.pdf", p5, width = 10, height = 8)
ggsave("figures/07_method_summary.png", p5, width = 10, height = 8, dpi = 150)

# 9. Top Hubs per Method (Bar Plot)

cat("Creating top hubs bar plot...\n")

top_10_per_method <- hub_metrics %>%
  group_by(method, gene) %>%
  summarise(total_degree = sum(degree), .groups = "drop") %>%
  group_by(method) %>%
  slice_max(total_degree, n = 10) %>%
  ungroup()

p6 <- ggplot(top_10_per_method, aes(x = reorder(gene, total_degree), 
                                     y = total_degree, fill = method)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~method, scales = "free_y") +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Top 10 Hub Genes per Method",
    subtitle = "Ranked by total degree across subtypes",
    x = "Gene",
    y = "Total Degree"
  ) +
  theme(legend.position = "none")

ggsave("figures/08_top_hubs_per_method.pdf", p6, width = 12, height = 8)
ggsave("figures/08_top_hubs_per_method.png", p6, width = 12, height = 8, dpi = 150)

# 10. Edge Overlap Visualization (UpSet-style)

cat("Creating edge overlap summary...\n")

# Simple bar plot of edges by method combination
edge_method_count <- read_csv("results/edge_method_count.csv", show_col_types = FALSE)

overlap_summary <- edge_method_count %>%
  count(methods, n_methods) %>%
  arrange(desc(n_methods), desc(n))

p7 <- ggplot(overlap_summary, aes(x = reorder(methods, n), y = n, fill = factor(n_methods))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  labs(
    title = "Edge Detection by Method Combination",
    subtitle = "Which methods detect which edges?",
    x = "Method Combination",
    y = "Number of Edges",
    fill = "# Methods"
  ) +
  theme(legend.position = "bottom")

ggsave("figures/09_edge_overlap.pdf", p7, width = 10, height = 6)
ggsave("figures/09_edge_overlap.png", p7, width = 10, height = 6, dpi = 150)

# 11. Summary Figure

cat("Creating summary figure...\n")

# Key metrics summary
p8 <- network_stats %>%
  ggplot(aes(x = n_nodes, y = n_edges, color = method, shape = subtype)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = method_colors) +
  labs(
    title = "Network Size: Nodes vs Edges",
    subtitle = "Each point = one method-subtype combination",
    x = "Number of Nodes",
    y = "Number of Edges",
    color = "Method",
    shape = "Subtype"
  ) +
  theme(legend.position = "right")

ggsave("figures/10_nodes_vs_edges.pdf", p8, width = 10, height = 6)
ggsave("figures/10_nodes_vs_edges.png", p8, width = 10, height = 6, dpi = 150)

# Done

cat("\n=== Visualization Complete ===\n")
cat("\nFigures saved in figures/ folder:\n")
cat("  01_network_size_comparison.pdf/png\n")
cat("  02_method_agreement_heatmap.pdf/png\n")
cat("  03_consensus_edges.pdf/png\n")
cat("  04_hub_gene_heatmap.pdf/png\n")
cat("  05_degree_distribution.pdf/png\n")
cat("  06_avg_degree_comparison.pdf/png\n")
cat("  07_method_summary.pdf/png\n")
cat("  08_top_hubs_per_method.pdf/png\n")
cat("  09_edge_overlap.pdf/png\n")
cat("  10_nodes_vs_edges.pdf/png\n")
