
# 06_visualize_causal.R
# Visualizations for causal edge analysis


library(tidyverse)
library(igraph)
library(ggplot2)
library(pheatmap)

cat("=== Causal Analysis Visualizations ===\n\n")

theme_set(theme_minimal(base_size = 12))

# Color palettes
subtype_colors <- c("basal" = "#E7298A", "her2" = "#7570B3", 
                    "luma" = "#1B9E77", "lumb" = "#D95F02")

# 1. Load Results

directed_edges <- read_csv("results/miic_directed_edges.csv", show_col_types = FALSE)
consensus_edges <- read_csv("results/consensus_3method_edges.csv", show_col_types = FALSE)
driver_edges <- read_csv("results/driver_gene_edges.csv", show_col_types = FALSE)
consensus_summary <- read_csv("results/consensus_summary.csv", show_col_types = FALSE)
regulators <- read_csv("results/miic_regulators.csv", show_col_types = FALSE)
targets <- read_csv("results/miic_targets.csv", show_col_types = FALSE)

cat("Loaded all causal analysis results\n")

# 2. Directed Edges per Subtype (Bar Plot)

cat("Creating directed edges bar plot...\n")

directed_summary <- directed_edges %>%
  group_by(subtype) %>%
  summarise(n_directed = n(), .groups = "drop")

p1 <- ggplot(directed_summary, aes(x = subtype, y = n_directed, fill = subtype)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n_directed), vjust = -0.5, size = 5) +
  scale_fill_manual(values = subtype_colors) +
  labs(
    title = "MIIC Directed Edges by Subtype",
    subtitle = "Edges with strong causal direction (p < 0.2)",
    x = "Subtype",
    y = "Number of Directed Edges"
  ) +
  theme(legend.position = "none") +
  ylim(0, max(directed_summary$n_directed) * 1.15)

ggsave("figures/causal_01_directed_edges.pdf", p1, width = 8, height = 6)
ggsave("figures/causal_01_directed_edges.png", p1, width = 8, height = 6, dpi = 150)

# 3. Three-Method Consensus (Stacked Bar)

cat("Creating consensus bar plot...\n")

consensus_long <- consensus_summary %>%
  select(subtype, two_methods, three_methods) %>%
  pivot_longer(-subtype, names_to = "consensus_level", values_to = "n_edges") %>%
  mutate(consensus_level = factor(consensus_level, 
                                   levels = c("two_methods", "three_methods"),
                                   labels = c("2 Methods", "3 Methods")))

p2 <- ggplot(consensus_long, aes(x = subtype, y = n_edges, fill = consensus_level)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("2 Methods" = "#6BAED6", "3 Methods" = "#08519C")) +
  labs(
    title = "Multi-Method Consensus Edges",
    subtitle = "High-confidence edges confirmed by multiple methods",
    x = "Subtype",
    y = "Number of Edges",
    fill = "Consensus Level"
  ) +
  theme(legend.position = "bottom")

ggsave("figures/causal_02_consensus_edges.pdf", p2, width = 8, height = 6)
ggsave("figures/causal_02_consensus_edges.png", p2, width = 8, height = 6, dpi = 150)

# 4. Top Regulatory Hubs Heatmap

cat("Creating regulatory hubs heatmap...\n")

# Get top 10 regulators per subtype
top_regs <- regulators %>%
  group_by(subtype) %>%
  slice_max(n_targets, n = 10) %>%
  ungroup()

# Pivot for heatmap
reg_wide <- top_regs %>%
  pivot_wider(names_from = subtype, values_from = n_targets, values_fill = 0) %>%
  column_to_rownames("source")

# Only plot if we have data
if (nrow(reg_wide) > 0 && ncol(reg_wide) > 0) {
  pdf("figures/causal_03_regulator_hubs.pdf", width = 8, height = 10)
  pheatmap(
    as.matrix(reg_wide),
    main = "Top Regulatory Hubs\n(Genes with most outgoing edges)",
    color = colorRampPalette(c("white", "#FEB24C", "#E31A1C"))(50),
    cluster_cols = FALSE,
    fontsize_row = 9,
    display_numbers = TRUE,
    number_format = "%d"
  )
  dev.off()
  
  png("figures/causal_03_regulator_hubs.png", width = 800, height = 1000)
  pheatmap(
    as.matrix(reg_wide),
    main = "Top Regulatory Hubs\n(Genes with most outgoing edges)",
    color = colorRampPalette(c("white", "#FEB24C", "#E31A1C"))(50),
    cluster_cols = FALSE,
    fontsize_row = 9,
    display_numbers = TRUE,
    number_format = "%d"
  )
  dev.off()
}

# 5. Top Target Hubs Heatmap

cat("Creating target hubs heatmap...\n")

# Get top 10 targets per subtype
top_targs <- targets %>%
  group_by(subtype) %>%
  slice_max(n_regulators, n = 10) %>%
  ungroup()

# Pivot for heatmap
targ_wide <- top_targs %>%
  pivot_wider(names_from = subtype, values_from = n_regulators, values_fill = 0) %>%
  column_to_rownames("target")

if (nrow(targ_wide) > 0 && ncol(targ_wide) > 0) {
  pdf("figures/causal_04_target_hubs.pdf", width = 8, height = 10)
  pheatmap(
    as.matrix(targ_wide),
    main = "Top Target Hubs\n(Genes with most incoming edges)",
    color = colorRampPalette(c("white", "#A1D99B", "#006D2C"))(50),
    cluster_cols = FALSE,
    fontsize_row = 9,
    display_numbers = TRUE,
    number_format = "%d"
  )
  dev.off()
  
  png("figures/causal_04_target_hubs.png", width = 800, height = 1000)
  pheatmap(
    as.matrix(targ_wide),
    main = "Top Target Hubs\n(Genes with most incoming edges)",
    color = colorRampPalette(c("white", "#A1D99B", "#006D2C"))(50),
    cluster_cols = FALSE,
    fontsize_row = 9,
    display_numbers = TRUE,
    number_format = "%d"
  )
  dev.off()
}

# 6. Driver Gene Network (per subtype)

cat("Creating driver gene networks...\n")

# Known drivers
known_drivers <- tolower(c(
  "ESR1", "GATA3", "FOXA1", "MYC", "E2F1", "E2F3", "E2F4", "RUNX1",
  "TP53", "PIK3CA", "ERBB2", "AKT1", "PTEN", "RB1", "CDK4", "CDK6", 
  "MAPK1", "MAPK3", "JAK1", "JAK2", "SRC", "CCND1", "CCNE1", "CCNB1", 
  "CDK2", "CDKN1A", "CDKN1B", "BCL2L1", "CASP8", "CASP9", "MMP1", "MMP11"
))

# Filter directed edges involving drivers
driver_directed <- directed_edges %>%
  filter(tolower(source) %in% known_drivers | tolower(target) %in% known_drivers)

if (nrow(driver_directed) > 0) {
  pdf("figures/causal_05_driver_networks.pdf", width = 12, height = 12)
  par(mfrow = c(2, 2))
  
  for (st in c("basal", "her2", "luma", "lumb")) {
    st_edges <- driver_directed %>% filter(subtype == st)
    
    if (nrow(st_edges) > 0) {
      g <- graph_from_data_frame(
        st_edges %>% select(source, target),
        directed = TRUE
      )
      
      # Color driver genes
      V(g)$color <- ifelse(tolower(V(g)$name) %in% known_drivers, 
                           subtype_colors[st], "lightgray")
      V(g)$label.cex <- 0.7
      
      plot(g, 
           main = paste(toupper(st), "- Driver Gene Causal Network"),
           vertex.size = 15,
           edge.arrow.size = 0.5,
           layout = layout_with_fr(g))
    }
  }
  dev.off()
}

# 7. Subtype-Specific Causal Edges Summary

cat("Creating subtype-specific summary...\n")

subtype_specific <- read_csv("results/subtype_specific_directed.csv", show_col_types = FALSE)

specific_counts <- subtype_specific %>%
  group_by(subtype) %>%
  summarise(n_specific = n(), .groups = "drop")

p3 <- ggplot(specific_counts, aes(x = subtype, y = n_specific, fill = subtype)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n_specific), vjust = -0.5, size = 5) +
  scale_fill_manual(values = subtype_colors) +
  labs(
    title = "Subtype-Specific Causal Edges",
    subtitle = "Directed edges unique to each subtype",
    x = "Subtype",
    y = "Number of Unique Directed Edges"
  ) +
  theme(legend.position = "none") +
  ylim(0, max(specific_counts$n_specific) * 1.15)

ggsave("figures/causal_06_subtype_specific.pdf", p3, width = 8, height = 6)
ggsave("figures/causal_06_subtype_specific.png", p3, width = 8, height = 6, dpi = 150)

# 8. Top Causal Edges Table (for report)

cat("Creating top causal edges summary table...\n")

# Top 5 strongest directed edges per subtype
top_causal <- directed_edges %>%
  group_by(subtype) %>%
  slice_max(direction_strength, n = 5) %>%
  ungroup() %>%
  select(subtype, source, target, direction_strength) %>%
  mutate(edge = paste(source, "→", target))

p4 <- ggplot(top_causal, aes(x = reorder(edge, direction_strength), 
                              y = direction_strength, fill = subtype)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~subtype, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = subtype_colors) +
  labs(
    title = "Top 5 Strongest Causal Edges per Subtype",
    subtitle = "Ranked by direction strength (|p_x2y - p_y2x|)",
    x = "",
    y = "Direction Strength"
  ) +
  theme(legend.position = "none")

ggsave("figures/causal_07_top_edges.pdf", p4, width = 12, height = 8)
ggsave("figures/causal_07_top_edges.png", p4, width = 12, height = 8, dpi = 150)

# 9. Summary Figure: Method Comparison

cat("Creating method comparison summary...\n")

method_comparison <- data.frame(
  Method = c("MIIC", "GENIE3", "ARACNE"),
  Type = c("Constraint-based\n(Causal)", "Tree-based\n(Predictive)", "MI-based\n(Association)"),
  Directed = c("Yes", "Yes", "No"),
  Total_Edges = c(sum(consensus_summary$miic), 
                  sum(consensus_summary$genie3), 
                  sum(consensus_summary$aracne))
)

p5 <- ggplot(method_comparison, aes(x = Method, y = Total_Edges, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(Total_Edges, "\n", Type)), 
            vjust = 1.2, color = "white", size = 3.5) +
  scale_fill_manual(values = c("MIIC" = "#E41A1C", "GENIE3" = "#377EB8", "ARACNE" = "#4DAF4A")) +
  labs(
    title = "Network Inference Methods Comparison",
    subtitle = "Total edges across all subtypes",
    x = "",
    y = "Total Edges"
  ) +
  theme(legend.position = "none")

ggsave("figures/causal_08_method_comparison.pdf", p5, width = 8, height = 6)
ggsave("figures/causal_08_method_comparison.png", p5, width = 8, height = 6, dpi = 150)

# 10. Combined Summary Figure

cat("Creating combined summary...\n")

# Key numbers for summary
summary_data <- data.frame(
  Metric = c("Total MIIC Edges", "Directed Edges", "3-Method Consensus", 
             "Driver Gene Edges", "Subtype-Specific"),
  Value = c(
    nrow(read_csv("results/miic_edges.csv", show_col_types = FALSE)),
    nrow(directed_edges),
    nrow(consensus_edges),
    nrow(driver_edges),
    nrow(subtype_specific)
  )
)

p6 <- ggplot(summary_data, aes(x = reorder(Metric, Value), y = Value)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_text(aes(label = Value), hjust = -0.2, size = 4) +
  coord_flip() +
  labs(
    title = "Causal Network Analysis Summary",
    x = "",
    y = "Number of Edges"
  ) +
  theme(axis.text.y = element_text(size = 11)) +
  xlim(NA, max(summary_data$Value) * 1.15)

ggsave("figures/causal_09_summary.pdf", p6, width = 10, height = 6)
ggsave("figures/causal_09_summary.png", p6, width = 10, height = 6, dpi = 150)

# Done

cat("\n=== Causal Visualizations Complete ===\n")
cat("\nFigures saved:\n")
cat("  causal_01_directed_edges.pdf/png\n")
cat("  causal_02_consensus_edges.pdf/png\n")
cat("  causal_03_regulator_hubs.pdf/png\n")
cat("  causal_04_target_hubs.pdf/png\n")
cat("  causal_05_driver_networks.pdf\n")
cat("  causal_06_subtype_specific.pdf/png\n")
cat("  causal_07_top_edges.pdf/png\n")
cat("  causal_08_method_comparison.pdf/png\n")
cat("  causal_09_summary.pdf/png\n")
