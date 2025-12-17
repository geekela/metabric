# Comparative Network Inference Across Breast Cancer Subtypes
## Benchmarking MIIC, GENIE3, and ARACNE on METABRIC

**MLBioNetw Course Project - CentraleSupélec**

## Biological Question

> Do different breast cancer molecular subtypes (Luminal A, Luminal B, HER2-enriched, Basal-like) exhibit distinct gene regulatory network architectures? Can we identify subtype-specific regulatory hubs that may serve as potential therapeutic targets?

## Project Overview

This project benchmarks three network inference methods on the METABRIC breast cancer gene expression dataset, stratified by PAM50 molecular subtype:

| Method | Type | Key Principle |
|--------|------|---------------|
| **MIIC** | Constraint-based | Conditional independence + information theory |
| **GENIE3** | Tree-based | Random forest variable importance |
| **ARACNE** | Information-theoretic | Mutual information + data processing inequality |

We compare the inferred networks across methods and subtypes to identify:
1. **Conserved vs. subtype-specific edges**
2. **Differential hub genes** with varying connectivity across subtypes
3. **Method agreement** on key regulatory interactions
4. **Biological validation** against known breast cancer driver genes

## Dataset

**METABRIC** (Molecular Taxonomy of Breast Cancer International Consortium)
- **Source**: [Kaggle](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric)
- **Samples**: ~1,904 breast cancer patients
- **Variables**: ~500 gene expression values (mRNA z-scores)
- **Subtypes**: Luminal A (~700), Luminal B (~400), HER2 (~200), Basal (~200)

## Project Structure

```
metabric/
├── data/                      # Raw and processed data (not tracked)
├── scripts/
│   ├── 01_preprocess.R        # Data loading and preprocessing
│   ├── 02a_run_miic.R         # MIIC network inference
│   ├── 02b_run_genie3.R       # GENIE3 network inference
│   ├── 02c_run_aracne.R       # ARACNE network inference
│   ├── 03_compare_networks.R  # Cross-method & cross-subtype comparison
│   └── 04_visualize.R         # Visualization and figures
├── results/                   # Output files (networks, metrics)
├── figures/                   # Generated plots
└── README.md
```

## Setup

### Required R Packages

```r
# From CRAN
install.packages(c("tidyverse", "miic", "igraph", "ggplot2", 
                   "pheatmap", "ppcor", "infotheo"))

# From Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GENIE3", "minet", "graph", "Rgraphviz"))

# Optional: for pathway enrichment
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
```

### Data Setup

1. Download from [Kaggle METABRIC](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric)
2. Place `METABRIC_RNA_Mutation.csv` in the `data/` folder

## Usage

```bash
# Run the full pipeline
Rscript scripts/01_preprocess.R        # ~1 min
Rscript scripts/02a_run_miic.R         # ~30-60 min
Rscript scripts/02b_run_genie3.R       # ~10-20 min
Rscript scripts/02c_run_aracne.R       # ~5-10 min
Rscript scripts/03_compare_networks.R  # ~2 min
Rscript scripts/04_visualize.R         # ~1 min
```

## Expected Outputs

### Results (`results/` folder)
- `miic_edges.csv`, `genie3_edges.csv`, `aracne_edges.csv` — Edge lists per method
- `all_edges_combined.csv` — Combined edge list with method labels
- `method_comparison.csv` — Cross-method agreement statistics
- `hub_comparison.csv` — Hub genes across methods and subtypes
- `driver_gene_analysis.csv` — Known driver genes in networks

### Figures (`figures/` folder)
- Network size comparison across methods/subtypes
- Jaccard similarity heatmaps (method × subtype)
- Hub gene agreement Venn diagrams
- Differential hub heatmaps
- Network visualizations

## Methods Summary

### MIIC (Multivariate Information-based Inductive Causation)
- Learns **causal/conditional independence** structure
- Can orient edges (infer causality direction)
- Handles latent confounders
- Reference: Verny et al. (2017) *PLOS Computational Biology*

### GENIE3 (GEne Network Inference with Ensemble of trees)
- **Random forest** regression for each gene
- Ranks regulators by variable importance
- Won DREAM4 network inference challenge
- Reference: Huynh-Thu et al. (2010) *PLOS ONE*

### ARACNE (Algorithm for the Reconstruction of Accurate Cellular Networks)
- **Mutual information** between gene pairs
- Applies **Data Processing Inequality** to remove indirect edges
- Widely used for regulatory network inference
- Reference: Margolin et al. (2006) *BMC Bioinformatics*

## References

1. Curtis C. et al. (2012). The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature*, 486(7403), 346-352.

2. Verny L. et al. (2017). Learning causal networks with latent variables from multivariate information in genomic data. *PLOS Computational Biology*, 13(10), e1005662.

3. Huynh-Thu V.A. et al. (2010). Inferring regulatory networks from expression data using tree-based methods. *PLOS ONE*, 5(9), e12776.

4. Margolin A.A. et al. (2006). ARACNE: An algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. *BMC Bioinformatics*, 7(Suppl 1), S7.

## Authors

MLBioNetw Course Project  
CentraleSupélec - AIDAMS Program  
Supervisors: Franck Simon, Hervé Isambert (Institut Curie)
