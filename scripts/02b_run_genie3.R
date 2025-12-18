
# 02b_run_genie3.R
# Run GENIE3 network inference on each breast cancer subtype


library(tidyverse)
library(GENIE3)

cat("=== GENIE3 Network Inference ===\n\n")


# 1. Load Preprocessed Data


expr_by_subtype <- readRDS("data/expr_by_subtype.rds")

cat("Loaded expression data for subtypes:\n")
for (st in names(expr_by_subtype)) {
  cat("  ", st, ":", nrow(expr_by_subtype[[st]]), "samples x", 
      ncol(expr_by_subtype[[st]]), "genes\n")
}


# 2. GENIE3 Parameters


# Configuration
N_TREES <- 100          # Number of trees per gene (default 1000, reduced for speed)
MAX_GENES <- 200        # Max genes to analyze
SAMPLE_CAP <- 500       # Max samples per subtype
N_CORES <- 4            # Parallel cores
TOP_EDGES <- 1000       # Keep top N edges per subtype

cat("\nParameters:\n")
cat("  Trees per gene:", N_TREES, "\n")
cat("  Max genes:", MAX_GENES, "\n")
cat("  Max samples:", SAMPLE_CAP, "\n")
cat("  Cores:", N_CORES, "\n")


# 3. Run GENIE3 on Each Subtype


run_genie3_subtype <- function(expr_df, subtype_name, n_trees = 100, 
                                max_genes = 200, sample_cap = 500, top_edges = 1000) {
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("Running GENIE3 on:", subtype_name, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  # Subsample if too many samples
  n_samples <- nrow(expr_df)
  if (n_samples > sample_cap) {
    set.seed(42)
    expr_df <- expr_df[sample(1:n_samples, sample_cap), ]
    cat("Subsampled to", sample_cap, "samples\n")
  }
  
  # Reduce genes if necessary
  n_genes <- ncol(expr_df)
  if (n_genes > max_genes) {
    gene_var <- apply(expr_df, 2, var)
    top_genes <- names(sort(gene_var, decreasing = TRUE))[1:max_genes]
    expr_df <- expr_df[, top_genes]
    cat("Reduced to top", max_genes, "genes by variance\n")
  }
  
  cat("Final dimensions:", nrow(expr_df), "samples x", ncol(expr_df), "genes\n")
  
  # GENIE3 expects genes as rows, samples as columns
  expr_matrix <- t(as.matrix(expr_df))
  
  cat("Running GENIE3 (this may take several minutes)...\n")
  start_time <- Sys.time()
  
  tryCatch({
    # Run GENIE3
    weight_matrix <- GENIE3(
      exprMatrix = expr_matrix,
      nCores = N_CORES,
      nTrees = n_trees,
      verbose = FALSE
    )
    
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "mins")
    cat("Completed in", round(elapsed, 2), "minutes\n")
    
    # Convert weight matrix to edge list
    edge_list <- getLinkList(weight_matrix, threshold = 0)
    colnames(edge_list) <- c("x", "y", "weight")
    
    # Keep top edges
    edge_list <- edge_list %>%
      arrange(desc(weight)) %>%
      head(top_edges)
    
    cat("Edges extracted:", nrow(edge_list), "\n")
    
    return(edge_list)
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(NULL)
  })
}

# Run on all subtypes
genie3_results <- list()

# Map to lowercase names to match MIIC online output
subtype_map <- c("LumA" = "luma", "LumB" = "lumb", "Her2" = "her2", "Basal" = "basal")
subtypes <- names(expr_by_subtype)

for (st in subtypes) {
  st_lower <- ifelse(st %in% names(subtype_map), subtype_map[st], tolower(st))
  genie3_results[[st_lower]] <- run_genie3_subtype(
    expr_by_subtype[[st]], 
    st_lower,
    n_trees = N_TREES,
    max_genes = MAX_GENES,
    sample_cap = SAMPLE_CAP,
    top_edges = TOP_EDGES
  )
}


# 4. Combine and Format Edges


cat("\n\n=== Edge Summary ===\n")

# Add subtype label and create edge_id
all_edges <- bind_rows(lapply(names(genie3_results), function(st) {
  edges <- genie3_results[[st]]
  if (is.null(edges) || nrow(edges) == 0) return(data.frame())
  
  edges %>%
    mutate(
      subtype = st,
      method = "GENIE3",
      # Create undirected edge ID for comparison
      edge_id = paste(pmin(x, y), pmax(x, y), sep = "---")
    )
}))

cat("\nTotal edges across all subtypes:", nrow(all_edges), "\n")
cat("Edges per subtype:\n")
print(table(all_edges$subtype))


# 5. Save Results


# Save raw results
saveRDS(genie3_results, "results/genie3_results.rds")

# Save edge table
write_csv(all_edges, "results/genie3_edges.csv")

# Save summary statistics
summary_stats <- data.frame(
  subtype = names(genie3_results),
  method = "GENIE3",
  n_edges = sapply(genie3_results, function(x) {
    if (is.null(x)) 0 else nrow(x)
  })
)
write_csv(summary_stats, "results/genie3_summary.csv")

cat("\n=== GENIE3 Inference Complete ===\n")
cat("Saved:\n")
cat("  - results/genie3_results.rds\n")
cat("  - results/genie3_edges.csv\n")
cat("  - results/genie3_summary.csv\n")
