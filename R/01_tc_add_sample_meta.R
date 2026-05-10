#' Add Sample-Level Metadata to a Seurat Object
#'
#' This function merges a data frame containing sample-level clinical or experimental
#' metadata into a Seurat object's cell-level metadata.
#'
#' @param object A Seurat object.
#' @param meta_df A data.frame containing the metadata.
#' @param join_by A character string specifying the column name to match.
#' @param factor_levels A named list specifying factor levels.
#' @param set_ident A character string to set as Idents.
#' @return A Seurat object with updated metadata.
#' @export
tc_add_sample_meta <- function(object, meta_df, join_by = "orig.ident", factor_levels = NULL, set_ident = NULL) {

  if (!inherits(meta_df, "data.frame")) stop("Error: meta_df must be a data.frame.")
  if (!join_by %in% colnames(meta_df)) stop(sprintf("Error: Column '%s' not found in meta_df.", join_by))
  if (!join_by %in% colnames(object@meta.data)) stop(sprintf("Error: Column '%s' not found in Seurat object.", join_by))

  curr_meta <- object@meta.data
  match_idx <- match(curr_meta[[join_by]], meta_df[[join_by]])

  if (any(is.na(match_idx))) {
    warning("Some cells have sample IDs that are not present in your meta_df.")
  }

  # 4. Inject new columns
  new_cols <- setdiff(colnames(meta_df), join_by)
  for (col in new_cols) {
    object@meta.data[[col]] <- meta_df[[col]][match_idx]
  }
  message(sprintf("✅ Successfully added %d metadata columns: %s", length(new_cols), paste(new_cols, collapse = ", ")))

  # 5. Apply factor levels
  if (!is.null(factor_levels) && is.list(factor_levels)) {
    for (col_name in names(factor_levels)) {
      if (col_name %in% colnames(object@meta.data)) {
        object@meta.data[[col_name]] <- factor(object@meta.data[[col_name]], levels = factor_levels[[col_name]])
        message(sprintf("✅ Converted '%s' to factor with levels: %s", col_name, paste(factor_levels[[col_name]], collapse = " -> ")))
      }
    }
  }

  # 6. Set Default Identity
  if (!is.null(set_ident)) {
    if (set_ident %in% colnames(object@meta.data)) {
      Seurat::Idents(object) <- set_ident
      message(sprintf("✅ Set active identity (Idents) to: '%s'", set_ident))
    }
  }

  return(object)
}
