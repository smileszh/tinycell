
# tinycell 🔬

**tinycell** is a lightweight, fully automated R wrapper package
specifically designed for 10X Genomics single-cell RNA sequencing
(scRNA-seq) data.

It condenses the lengthy standard Seurat and scDblFinder workflows into
highly streamlined `tc_` prefix functions. Its goal is to help
researchers complete single-cell data pre-processing, doublet removal,
and batch-effect-corrected dimensionality reduction with minimal code
and maximum speed.

## 🌟 Core Features

- **Seamless Import**: Automatically scan and merge multiple 10X sample
  directories with a single command (`tc_read_10x`).
- **Smart Quality Control**: Integrated bidirectional filtering based on
  mitochondrial content, features (genes), and UMIs (`tc_qc_filter`).
- **Ultra-Fast Doublet Removal**: Natively powered by the
  next-generation `scDblFinder` algorithm, featuring automated
  multi-sample parallelization and Seurat v5 Assay layer compatibility
  (`tc_remove_doublets`).
- **One-Click Integration**: Incorporates the `Harmony` algorithm to
  perform the entire pipeline—from normalization to UMAP—in one line of
  code (`tc_process_harmony`).

## 📥 Installation

You can install the development version from GitHub:

``` r
# Install devtools if you haven't already: install.packages("devtools")
devtools::install_github("smileszh/tinycell")
```

## 🚀 Quick Start

Experience the joy of compressing an 80-line pipeline script into just 4
lines of elegant code:

``` r
library(tinycell)

# 1. Automatically read and merge all 10X samples in a directory
sce <- tc_read_10x(data_dirs = "./data/", project_name = "My_scRNA_Project")

# 2. Quality Control & Filtering (e.g., genes > 300, MT < 20%)
sce <- tc_qc_filter(sce, min_features = 300, max_mt = 20)

# 3. Automatically identify and remove doublets (processed independently by sample)
sce <- tc_remove_doublets(sce, samples = "orig.ident")

# 4. Fully automated dimensionality reduction & Harmony batch integration
sce <- tc_process_harmony(sce, batch_var = "orig.ident")

# Inspect your beautifully cleaned and integrated object!
print(sce)
```

## 🤝 Dependencies & Acknowledgements

`tinycell` stands on the shoulders of giants. We deeply rely on and are
grateful for the following outstanding open-source projects: \*
[Seurat](https://satijalab.org/seurat/) (Core data infrastructure) \*
[scDblFinder](https://bioconductor.org/packages/scDblFinder/) (Doublet
identification) \* [Harmony](https://github.com/immunogenomics/harmony)
(Batch effect integration)
