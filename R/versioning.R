# ==============================================================================
# versioning.R — Auto-versioned output folder management
# ==============================================================================

#' Create an auto-versioned folder (e.g., GSA_v1, GSA_v2, ...)
#'
#' @param base_path Base directory where versioned folders are created
#' @param prefix    Version prefix (e.g., "GSA_v")
#' @return List with $tag (e.g., "GSA_v3") and $path (full path)
create_versioned_folder <- function(base_path, prefix) {
  folders  <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)
  versions <- suppressWarnings(
    as.numeric(gsub(prefix, "", folders[grep(paste0("^", prefix), folders)]))
  )

  next_v <- if (length(versions) == 0 || all(is.na(versions))) 1L else max(versions, na.rm = TRUE) + 1L

  tag  <- paste0(prefix, next_v)
  path <- file.path(base_path, tag)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  list(tag = tag, path = path)
}

#' Find the latest versioned folder matching a prefix
#'
#' @param base_path Base directory to search in
#' @param prefix    Folder prefix (e.g., "GSA_v")
#' @return List with $tag and $path of the latest version
find_latest_folder <- function(base_path, prefix) {
  folders  <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)
  matches  <- grep(paste0("^", prefix), folders, value = TRUE)

  if (length(matches) == 0) {
    stop("No folders found with prefix '", prefix, "' in: ", base_path)
  }

  versions <- suppressWarnings(as.numeric(gsub(prefix, "", matches)))
  latest_v <- max(versions, na.rm = TRUE)
  tag      <- paste0(prefix, latest_v)

  list(tag = tag, path = file.path(base_path, tag))
}
