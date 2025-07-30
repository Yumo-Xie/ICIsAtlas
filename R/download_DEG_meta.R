#' @importFrom rappdirs           user_data_dir

NULL

#' Download and load the DEG_meta dataset
#'
#' This function checks for a cached copy of DEG_meta.rda, downloads it from the
#' GitHub release if missing, and returns the DEG_meta object.
#'
#' @param destfile Optional path to store the downloaded file. Defaults to a
#'   user-specific cache directory via rappdirs.
#' @return The DEG_meta object as a data.frame or tibble (whatever is inside the .rda).
#' @export
download_DEG_meta <- function(destfile = NULL) {
  url <- "https://github.com/Yumo-Xie/ICIsAtlas/releases/download/v0.1.0/DEG_meta.rda"

  if (exists("DEG_meta", envir = .GlobalEnv)) {
    # message("DEG_meta already loaded in memory.")
    return(invisible(NULL))
  }

  # Decide where to store the file
  if (is.null(destfile)) {
    destfile <- file.path(rappdirs::user_data_dir("ICIsAtlas"), "DEG_meta.rda")
  }
  if (!dir.exists(dirname(destfile))) dir.create(dirname(destfile), recursive = TRUE)

  # Download if not cached
  if (!file.exists(destfile)) {
    message("Downloading DEG_meta.rda (~327 MB)...")
    utils::download.file(url, destfile, mode = "wb")
  } else {
    message("Using cached DEG_meta.rda")
  }

  # Load and return the object
  e <- new.env()
  load(destfile, envir = e)
  assign("DEG_meta", e$DEG_meta, envir = .GlobalEnv)
  invisible(NULL)
}
