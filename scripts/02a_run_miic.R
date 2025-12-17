# =============================================================================
# 02_run_miic.R
# Run MIIC network inference on each breast cancer subtype
# =============================================================================

library(tidyverse)
library(miic)

cat("=== MIIC Network Inference ===\n\n")

# -----------------------------------------------------------------------------
# 1. Load Preprocessed Data
# -----------------------------------------------------------------------------

expr_by_subtype <- readRDS("data/expr_by_subtype.rds")

cat("Loaded expression data for subtypes:\n")
for (st in names(expr_by_subtype)) {
  cat("  ", st, ":", nrow(expr_by_subtype[[st]]), "samples x", 
      ncol(expr_by_subtype[[st]]), "genes\n")
}

# -----------------------------------------------------------------------------
# 2. MIIC Parameters
# -----------------------------------------------------------------------------

# Configuration
N_THREADS <- 4          # Adjust based on your CPU
MAX_GENES <- 200        # Reduce if MIIC is too slow
SAMPLE_CAP <- 500       # Max samples per subtype (speeds up computation)

cat("\nParameters:\n")
cat("  Max genes:", MAX_GENES, "\n")
cat("  Max samples per subtype:", SAMPLE_CAP, "\n")
cat("  Threads:", N_THREADS, "\n")

# -----------------------------------------------------------------------------
# 3. Run MIIC on Each Subtype
# -----------------------------------------------------------------------------

run_miic_subtype <- function(expr_df, subtype_name, max_genes = 200, sample_cap = 500) {
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("Running MIIC on:", subtype_name, "\n")
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
  
  cat("Final dimensions:", nrow(expr_df), "x", ncol(expr_df), "\n")
  
  # Prepare data for MIIC
  data_miic <- as.data.frame(expr_df)
  
  # Run MIIC
  cat("Running MIIC (this may take several minutes)...\n")
  start_time <- Sys.time()
  
  tryCatch({
    result <- miic(
      input_data = data_miic,
      cplx = "nml",           # Normalized Maximum Likelihood complexity
      ori = TRUE,             # Enable edge orientation
      propagation = FALSE,    # Don't propagate orientations
      n_threads = N_THREADS
    )
    
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "mins")
    cat("Completed in", round(elapsed, 2), "minutes\n")
    
    # Summary
    if (!is.null(result$all.edges.summary)) {
      n_edges <- nrow(result$all.edges.summary)
      cat("Edges found:", n_edges, "\n")
    } else {
      cat("Warning: No edges found\n")
    }
    
    return(result)
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(NULL)
  })
}

# Run on all subtypes
miic_results <- list()
subtypes <- names(expr_by_subtype)

for (st in subtypes) {
  miic_results[[st]] <- run_miic_subtype(
    expr_by_subtype[[st]], 
    st,
    max_genes = MAX_GENES,
    sample_cap = SAMPLE_CAP
  )
}

# -----------------------------------------------------------------------------
# 4. Extract Edge Summaries
# -----------------------------------------------------------------------------

cat("\n\n=== Edge Summary ===\n")

extract_edges <- function(miic_result, subtype_name) {
  if (is.null(miic_result) || is.null(miic_result$all.edges.summary)) {
    return(data.frame())
  }
  
  edges <- miic_result$all.edges.summary
  
  edges %>%
    mutate(
      subtype = subtype_name,
      method = "MIIC",
      # Create undirected edge ID for comparison
      edge_id = paste(pmin(x, y), pmax(x, y), sep = "---")
    )
}

all_edges <- bind_rows(lapply(names(miic_results), function(st) {
  extract_edges(miic_results[[st]], st)
}))

cat("\nTotal edges across all subtypes:", nrow(all_edges), "\n")
cat("Edges per subtype:\n")
print(table(all_edges$subtype))

# -----------------------------------------------------------------------------
# 5. Save Results
# -----------------------------------------------------------------------------

# Save MIIC results
saveRDS(miic_results, "results/miic_results.rds")

# Save edge table
write_csv(all_edges, "results/miic_edges.csv")

# Save summary statistics
summary_stats <- data.frame(
  subtype = names(miic_results),
  method = "MIIC",
  n_edges = sapply(names(miic_results), function(st) {
    if (is.null(miic_results[[st]]$all.edges.summary)) 0 
    else nrow(miic_results[[st]]$all.edges.summary)
  }),
  n_samples = sapply(names(expr_by_subtype), function(st) {
    min(nrow(expr_by_subtype[[st]]), SAMPLE_CAP)
  })
)
write_csv(summary_stats, "results/miic_summary.csv")

cat("\n=== MIIC Inference Complete ===\n")
cat("Saved:\n")
cat("  - results/miic_results.rds\n")
cat("  - results/miic_edges.csv\n")
cat("  - results/miic_summary.csv\n")
