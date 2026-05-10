#' Quality Control and Filtering of a Seurat Object
#'
#' This function calculates mitochondrial gene percentage and filters cells based on
#' features (genes), total counts (UMIs), and mitochondrial content.
#'
#' @param object A Seurat object.
#' @param mt_pattern Regex pattern for mitochondrial genes. Default is \code{"^MT-"}.
#' @param min_features Min genes per cell. Default is 300.
#' @param max_features Max genes per cell. Default is 7500.
#' @param min_counts Min total counts (UMI) per cell. Default is 500.
#' @param max_counts Max total counts (UMI) per cell. Default is Inf (no upper limit).
#' @param max_mt Max percentage of mitochondrial genes. Default is 20.
#'
#' @return A filtered \code{Seurat} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sce.filt <- tc_qc_filter(
#'   object = sce.all,
#'   min_features = 300,
#'   max_features = 7500,
#'   min_counts = 500,
#'   max_mt = 20
#' )
#' }
tc_qc_filter <- function(object,
                         mt_pattern = "^MT-",
                         min_features = 300,
                         max_features = 7500,
                         min_counts = 500,
                         max_counts = Inf,
                         max_mt = 20) {

  # 1. Calculate percent.mt
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = mt_pattern)

  n_cells_before <- ncol(object)
  meta <- object@meta.data

  # 2. Refined Filtering Logic
  # Added nCount_RNA conditions
  cells_to_keep <- rownames(meta)[
    meta$nFeature_RNA > min_features &
      meta$nFeature_RNA < max_features &
      meta$nCount_RNA   > min_counts   &
      meta$nCount_RNA   < max_counts   &
      meta$percent.mt   < max_mt
  ]

  object <- subset(object, cells = cells_to_keep)

  # 3. Stats reporting
  n_cells_after <- ncol(object)
  n_removed <- n_cells_before - n_cells_after

  message("✅ QC Filtering Summary:")
  message(sprintf("   ➤ Original cells : %d", n_cells_before))
  message(sprintf("   ➤ Filtered cells : %d", n_cells_after))
  message(sprintf("   ➤ Removed cells  : %d (%.2f%%)", n_removed, (n_removed/n_cells_before)*100))

  return(object)
}
