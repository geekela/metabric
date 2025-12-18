
# 02c_run_aracne.R
# Run ARACNE network inference on each breast cancer subtype
# Uses the minet package implementation


library(tidyverse)
library(minet)
library(infotheo)

cat("=== ARACNE Network Inference ===\n\n")


# 1. Load Preprocessed Data


expr_by_subtype <- readRDS("data/expr_by_subtype.rds")

cat("Loaded expression data for subtypes:\n")
for (st in names(expr_by_subtype)) {
  cat("  ", st, ":", nrow(expr_by_subtype[[st]]), "samples x", 
      ncol(expr_by_subtype[[st]]), "genes\n")
}

# 2. ARACNE Parameters  

# Configuration
MAX_GENES <- 200        # Max genes to analyze
SAMPLE_CAP <- 500       # Max samples per subtype
N_BINS <- 10            # Discretization bins for MI estimation
TOP_EDGES <- 1000       # Keep top N edges per subtype
MI_ESTIMATOR <- "mi.shrink"  # MI estimator: "mi.empirical", "mi.mm", "mi.shrink"

cat("\nParameters:\n")
cat("  Max genes:", MAX_GENES, "\n")
cat("  Max samples:", SAMPLE_CAP, "\n")
cat("  Discretization bins:", N_BINS, "\n")
cat("  MI estimator:", MI_ESTIMATOR, "\n")

# 3. Run ARACNE on Each Subtype

run_aracne_subtype <- function(expr_df, subtype_name, max_genes = 200, 
                                sample_cap = 500, n_bins = 10, top_edges = 1000) {
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("Running ARACNE on:", subtype_name, "\n")
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
  
  # minet expects samples as rows, genes as columns
  expr_matrix <- as.matrix(expr_df)
  
  # Discretize for mutual information estimation
  cat("Discretizing expression data...\n")
  expr_disc <- discretize(expr_matrix, disc = "equalfreq", nbins = n_bins)
  
  cat("Running ARACNE (this may take a few minutes)...\n")
  start_time <- Sys.time()
  
  tryCatch({
    # Build mutual information matrix
    mi_matrix <- build.mim(expr_disc, estimator = MI_ESTIMATOR)
    
    # Apply ARACNE (Data Processing Inequality)
    aracne_matrix <- aracne(mi_matrix, eps = 0.05)
    
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "mins")
    cat("Completed in", round(elapsed, 2), "minutes\n")
    
    # Convert adjacency matrix to edge list
    genes <- colnames(aracne_matrix)
    edge_list <- data.frame()
    
    for (i in 1:(length(genes) - 1)) {
      for (j in (i + 1):length(genes)) {
        weight <- aracne_matrix[i, j]
        if (weight > 0) {
          edge_list <- rbind(edge_list, data.frame(
            x = genes[i],
            y = genes[j],
            weight = weight
          ))
        }
      }
    }
    
    # Sort by weight and keep top edges
    if (nrow(edge_list) > 0) {
      edge_list <- edge_list %>%
        arrange(desc(weight)) %>%
        head(top_edges)
    }
    
    cat("Edges extracted:", nrow(edge_list), "\n")
    
    return(edge_list)
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(NULL)
  })
}

# Run on all subtypes
aracne_results <- list()

# Map to lowercase names to match MIIC online output
subtype_map <- c("LumA" = "luma", "LumB" = "lumb", "Her2" = "her2", "Basal" = "basal")
subtypes <- names(expr_by_subtype)

for (st in subtypes) {
  st_lower <- ifelse(st %in% names(subtype_map), subtype_map[st], tolower(st))
  aracne_results[[st_lower]] <- run_aracne_subtype(
    expr_by_subtype[[st]], 
    st_lower,
    max_genes = MAX_GENES,
    sample_cap = SAMPLE_CAP,
    n_bins = N_BINS,
    top_edges = TOP_EDGES
  )
}

# 4. Combine and Format Edges 

cat("\n\n=== Edge Summary ===\n")

# Add subtype label and create edge_id
all_edges <- bind_rows(lapply(names(aracne_results), function(st) {
  edges <- aracne_results[[st]]
  if (is.null(edges) || nrow(edges) == 0) return(data.frame())

  edges %>%
    mutate(
      subtype = st,
      method = "ARACNE",
      # Create undirected edge ID for comparison
      edge_id = paste(pmin(x, y), pmax(x, y), sep = "---")
    )
}))

cat("\nTotal edges across all subtypes:", nrow(all_edges), "\n")
cat("Edges per subtype:\n")
print(table(all_edges$subtype))


# 5. Save Results


# Save raw results
saveRDS(aracne_results, "results/aracne_results.rds")

# Save edge table
write_csv(all_edges, "results/aracne_edges.csv")

# Save summary statistics
summary_stats <- data.frame(
  subtype = names(aracne_results),
  method = "ARACNE",
  n_edges = sapply(aracne_results, function(x) {
    if (is.null(x)) 0 else nrow(x)
  })
)
write_csv(summary_stats, "results/aracne_summary.csv")

cat("\n=== ARACNE Inference Complete ===\n")
cat("Saved:\n")
cat("  - results/aracne_results.rds\n")
cat("  - results/aracne_edges.csv\n")
cat("  - results/aracne_summary.csv\n")
