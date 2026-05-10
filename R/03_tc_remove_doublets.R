#' Identify and Remove Doublets using scDblFinder
#'
#' This function uses the fast and modern \code{scDblFinder} algorithm.
#' It includes native compatibility with Seurat v5 Assay layers, safely
#' joining them in a temporary environment without disrupting the downstream
#' v5 integration workflow.
#'
#' @param object A Seurat object.
#' @param samples Metadata column containing sample identities. Default is \code{"orig.ident"}.
#' @param filter_cells Logical. If TRUE, removes doublets. Default is \code{TRUE}.
#'
#' @return A Seurat object cleaned of doublets.
#' @export
#'
#' @importFrom SingleCellExperiment colData
#'
tc_remove_doublets <- function(object, samples = "orig.ident", filter_cells = TRUE) {

  if (!requireNamespace("scDblFinder", quietly = TRUE)) stop("Install scDblFinder via BiocManager")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("Install SingleCellExperiment via BiocManager")
  if (!samples %in% colnames(object@meta.data)) stop(sprintf("Column '%s' not found.", samples))

  message("🚀 Preparing Seurat object for conversion...")

  # ==========================================
  # Seurat v5 Compatibility Layer
  # ==========================================
  temp_obj <- object
  # Check if it's a v5 assay and has split layers
  if (inherits(temp_obj[["RNA"]], "Assay5")) {
    message("   ➤ Detected Seurat v5 Assay. Temporarily joining layers for conversion...")
    temp_obj <- SeuratObject::JoinLayers(temp_obj)
  }

  message("🚀 Converting to SingleCellExperiment...")
  sce <- suppressWarnings(Seurat::as.SingleCellExperiment(temp_obj))

  # Free up memory immediately
  rm(temp_obj)
  gc()

  # ==========================================
  # Run scDblFinder
  # ==========================================
  message(sprintf("🧬 Running scDblFinder (auto-processing by '%s')...", samples))
  sce <- scDblFinder::scDblFinder(sce, samples = samples)

  # ==========================================
  # Transfer Results back to ORIGINAL object
  # ==========================================
  message("🔄 Transferring results back to Seurat object...")
  dbl_class <- SingleCellExperiment::colData(sce)$scDblFinder.class
  dbl_score <- SingleCellExperiment::colData(sce)$scDblFinder.score

  # Inject metadata into the original object (preserving v5 split layers!)
  object$doublet_finder <- ifelse(dbl_class == "doublet", "Doublet", "Singlet")
  object$doublet_score  <- dbl_score

  dbl_table <- table(object$doublet_finder)
  message("✅ Doublet detection complete!")
  print(dbl_table)

  if (filter_cells) {
    message("🧹 Removing doublets...")
    cells_to_keep <- rownames(object@meta.data)[object@meta.data$doublet_finder == "Singlet"]
    object <- subset(object, cells = cells_to_keep)
    message(sprintf("   Cells remaining: %d", ncol(object)))
  }

  return(object)
}
