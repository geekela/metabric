
# 02a_load_miic.R
# Load MIIC results from online platform (https://miic.curie.fr)

library(tidyverse)

cat("=== Loading MIIC Online Results ===\n\n")


# 1. Define Subtypes and Paths


subtypes <- c("basal", "luma", "lumb", "her2")
miic_dir <- "results/miic_online"

# Check that directories exist
for (st in subtypes) {
  path <- file.path(miic_dir, st)
  if (!dir.exists(path)) {
    stop(paste("Directory not found:", path))
  }
}

cat("Found MIIC results for subtypes:", paste(subtypes, collapse = ", "), "\n")


# 2. Load Edge Lists from Each Subtype


load_miic_edges <- function(subtype) {
  edge_file <- file.path(miic_dir, subtype, "edgesList.miic.summary.txt")
  
  if (!file.exists(edge_file)) {
    warning(paste("Edge file not found:", edge_file))
    return(data.frame())
  }
  
  cat("Loading:", edge_file, "\n")
  
  edges <- read.delim(edge_file, stringsAsFactors = FALSE)
  
  cat("  Found", nrow(edges), "edges\n")
  
  # Standardize and add metadata
  edges %>%
    mutate(
      subtype = subtype,
      method = "MIIC",
      # Create undirected edge ID for comparison
      edge_id = paste(pmin(x, y), pmax(x, y), sep = "---")
    ) %>%
    select(x, y, edge_id, subtype, method, everything())
}

# Load all subtypes
miic_edges_list <- lapply(subtypes, load_miic_edges)
names(miic_edges_list) <- subtypes

# Combine into single dataframe
all_edges <- bind_rows(miic_edges_list)


# 3. Summary Statistics


cat("\n=== MIIC Edge Summary ===\n")
cat("Total edges across all subtypes:", nrow(all_edges), "\n\n")

cat("Edges per subtype:\n")
print(table(all_edges$subtype))

# Unique genes
genes <- unique(c(all_edges$x, all_edges$y))
cat("\nTotal unique genes:", length(genes), "\n")

# Edge overlap preview
edge_counts <- all_edges %>%
  group_by(edge_id) %>%
  summarise(n_subtypes = n_distinct(subtype), .groups = "drop")

cat("\nEdge conservation:\n")
print(table(edge_counts$n_subtypes))


# 4. Save Standardized Results


# Save combined edge table
write_csv(all_edges, "results/miic_edges.csv")

# Save summary
summary_stats <- data.frame(
  subtype = subtypes,
  method = "MIIC",
  n_edges = sapply(miic_edges_list, nrow)
)
write_csv(summary_stats, "results/miic_summary.csv")

cat("\n=== MIIC Loading Complete ===\n")
cat("Saved:\n")
cat("  - results/miic_edges.csv\n")
cat("  - results/miic_summary.csv\n")
