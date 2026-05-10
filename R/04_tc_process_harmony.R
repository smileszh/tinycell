#' Run Standard Integration and Dimensionality Reduction via Harmony
#'
#' This function encapsulates the standard single-cell processing pipeline:
#' Normalization, Feature Selection, Scaling, PCA, Harmony batch correction,
#' Graph-based Clustering, and non-linear dimensionality reduction (UMAP/tSNE).
#'
#' @param object A Seurat object.
#' @param batch_var Character. The metadata column to use for Harmony integration
#'   (batch effect removal). Default is \code{"orig.ident"}.
#' @param dims Numeric vector. The PCA dimensions to use. Default is \code{1:30}.
#' @param resolution Numeric or numeric vector. Resolution parameter(s) for clustering.
#'   If a vector is provided, it calculates all resolutions (useful for clustree).
#'   Default is \code{0.3}.
#' @param run_tsne Logical. Whether to run tSNE alongside UMAP. Default is \code{FALSE}
#'   to save computation time.
#'
#' @return A Seurat object with PCA, Harmony, UMAP (and optionally tSNE) reductions,
#'   along with clustering results.
#' @export
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters RunUMAP RunTSNE
tc_process_harmony <- function(object, batch_var = "orig.ident", dims = 1:30, resolution = 0.3, run_tsne = FALSE) {

  if (!requireNamespace("harmony", quietly = TRUE)) stop("Please install the 'harmony' package.")
  if (!batch_var %in% colnames(object@meta.data)) stop(sprintf("Column '%s' not found.", batch_var))

  message("🧪 Step 1: Normalization, Feature Selection, and Scaling...")
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, verbose = FALSE)
  object <- Seurat::ScaleData(object, verbose = FALSE)

  message("📉 Step 2: Running PCA...")
  object <- Seurat::RunPCA(object, verbose = FALSE)

  message(sprintf("🧬 Step 3: Running Harmony integration grouped by '%s'...", batch_var))
  object <- harmony::RunHarmony(object, group.by.vars = batch_var, verbose = FALSE)

  message(sprintf("🕸️ Step 4: Building Graph using Harmony dimensions %d to %d...", min(dims), max(dims)))
  object <- Seurat::FindNeighbors(object, reduction = "harmony", dims = dims, verbose = FALSE)

  message(sprintf("📊 Step 5: Finding Clusters at resolution(s): %s...", paste(resolution, collapse = ", ")))
  object <- Seurat::FindClusters(object, resolution = resolution, verbose = FALSE)

  message("🗺️  Step 6: Running UMAP...")
  object <- Seurat::RunUMAP(object, reduction = "harmony", dims = dims, verbose = FALSE)

  if (run_tsne) {
    message("🗺️  Step 7: Running tSNE...")
    object <- Seurat::RunTSNE(object, reduction = "harmony", dims = dims, verbose = FALSE)
  }

  message("✅ Integration and Dimensionality Reduction Complete!")
  return(object)
}
