#' Read and Merge Single or Multiple 10X Genomics Datasets
#'
#' This function provides a streamlined workflow for reading 10X Genomics data and
#' creating a merged Seurat object. It supports two input modes:
#' 1. A single parent directory path: The function will automatically detect all
#'    subdirectories containing 10X data and merge them.
#' 2. A named character vector: The function will read the specific paths provided
#'    and use the vector names as sample identifiers.
#'
#' @param data_dirs A single character string specifying the parent directory, OR a
#'   named character vector of specific sample paths (e.g., \code{c(Ctrl = "./path1")}).
#' @param project_name Character string to set the project name of the merged Seurat
#'   object. Default is \code{"tinycell_project"}.
#' @param min_cells Include features detected in at least this many cells. Default is 3.
#' @param min_features Include cells where at least this many features are detected. Default is 200.
#'
#' @return A merged \code{Seurat} object containing all processed samples. Cell names
#'   are automatically prefixed with their respective sample names to prevent barcode collisions.
#' @export
#'
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' # Mode 1: Auto-detect all sample folders inside a parent directory
#' my_seurat <- tc_read_10x(data_dirs = "./data/raw_samples/")
#'
#' # Mode 2: Specify paths manually using a named vector
#' paths <- c(Sample_A = "./data/SRR1", Sample_B = "./data/SRR2")
#' my_seurat <- tc_read_10x(data_dirs = paths, project_name = "My_Experiment")
#' }
tc_read_10x <- function(data_dirs, project_name = "tinycell_project", min_cells = 3, min_features = 200) {

  # ==========================================
  # Smart Detection: If input is a single parent directory
  # ==========================================
  if (length(data_dirs) == 1 && is.null(names(data_dirs))) {
    if (!dir.exists(data_dirs)) {
      stop(sprintf("Error: Directory '%s' does not exist.", data_dirs))
    }

    # Auto-detect all subdirectories
    sub_dirs <- list.dirs(data_dirs, full.names = TRUE, recursive = FALSE)

    if (length(sub_dirs) == 0) {
      # The directory itself might be a single 10X output folder
      data_dirs <- setNames(data_dirs, basename(data_dirs))
    } else {
      # Convert to a named vector using folder names
      data_dirs <- setNames(sub_dirs, basename(sub_dirs))
      message(sprintf("🤖 Auto-detected %d sample directories. Starting batch processing...", length(data_dirs)))
    }
  }

  # ==========================================
  # Core Reading and Object Creation
  # ==========================================
  seurat_list <- list()

  for (sample_name in names(data_dirs)) {
    path <- data_dirs[[sample_name]]
    message(sprintf("➤ Reading sample: %s ...", sample_name))

    # Read 10X matrix
    counts <- Seurat::Read10X(data.dir = path)

    # Print raw dimensions
    message(sprintf("   Matrix dimensions: %d features x %d cells", nrow(counts), ncol(counts)))

    # Create Seurat object
    obj <- Seurat::CreateSeuratObject(counts = counts,
                                      project = sample_name,
                                      min.cells = min_cells,
                                      min.features = min_features)
    seurat_list[[sample_name]] <- obj
  }

  # ==========================================
  # Merging Process
  # ==========================================
  if (length(seurat_list) == 1) {
    message("✅ Processing complete (Single Sample).")
    return(seurat_list[[1]])
  } else {
    message("🔄 Merging all samples...")
    base_obj <- seurat_list[[1]]
    rest_objs <- seurat_list[-1]

    # Merge and prepend cell IDs with sample names to avoid barcode collisions
    merged_obj <- merge(x = base_obj,
                        y = rest_objs,
                        add.cell.ids = names(seurat_list),
                        project = project_name)

    message("✅ All samples merged successfully!")
    return(merged_obj)
  }
}
