#' Recover cell‐states and FavTME from bulk expression
#'
#' Runs the EcoTyper pipeline on a expression matrix, returning both
#' cell-state abundances and a composite tumor-microenvironment favorability score.
#'
#' @param expr Numeric matrix or data.frame; genes in rows and samples in columns.
#' @param expr_type Character; one of \code{"array"}, \code{"counts"}, \code{"FPKM"}, or \code{"TPM"} (default: \code{"array"}).
#' @param threads Integer; number of parallel threads to use (default: \code{1}).
#' @param outdir Character; path to the output directory (default: \code{tempdir()}).
#' @param ecotyper_dir Character; path to the directory containing EcoTyper scripts.
#' @return A \code{list} with components \code{FavTME} (favorability score) and
#'   \code{CellStates} (cell-state abundance table).
#' @examples
#' \dontrun{
#' data(example_mat, package = "ICIsAtlas")
#' res <- recover_cellstates(
#'   expr        = example_mat,
#'   expr_type   = "counts",
#'   threads     = 2
#' )
#' head(res$FavTME)
#' }
#' @export
recover_cellstates <- function(expr,
                             expr_type    = c("array", "counts", "FPKM", "TPM"),
                             threads      = 1,
                             outdir       = tempdir(),
                             ecotyper_dir = system.file("extdata/ecotyper", package = "ICIsAtlas")
                             ) {
  expr_type <- match.arg(expr_type)


  if(expr_type == "counts"){
    expr = cal_rpkm(expr)
  }

  expr = cbind(rownames(expr),expr)
  colnames(expr)[1] = "Gene"

  if (is.matrix(expr) || is.data.frame(expr)) {
    expr_file <- file.path(outdir, "expression_input.txt")
    write.table(expr, expr_file, sep = "\t", quote = FALSE, row.names  = F)
  } else {
    stop("`expr` must be a matrix/data.frame.")
  }

  if (.Platform$OS.type == "windows") {
    rscript <- normalizePath(file.path(R.home("bin"), "Rscript.exe"),
                             winslash = "/", mustWork = TRUE)
  } else {
    rscript <- Sys.which("Rscript")
    if (!nzchar(rscript)) stop("`Rscript` not found on PATH.")
  }


  script_name <- "EcoTyper_recovery_bulk.R"
  args <- c(
    script_name,
    "-d", "Carcinoma",
    "-m", normalizePath(expr_file, winslash = "/", mustWork = TRUE),
    "-o", normalizePath(outdir, winslash = "/", mustWork = TRUE),
    "-t", threads
  )

    args <- shQuote(args, type = if (.Platform$OS.type=="windows") "cmd" else "sh")


  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(ecotyper_dir)

  log <- system2(rscript, args = args, stdout = TRUE, stderr = TRUE)
  cat(paste(log, collapse = "\n"), "\n")

  status <- attr(log, "status")
  if (!is.null(status) && status != 0) {
    stop("EcoTyper pipeline failed (exit code ", status, "). See logs above.")
  }


  key_dir = file.path("./EcoTyper", "Carcinoma", "Carcinoma_Fractions", "Analysis", "rank_selection")
  states_dir = file.path(outdir,"expression_input")
  key = read.delim(file.path(key_dir, "rank_data.txt"))
  inpput_dir = file.path("../EcoTyper", "Carcinoma", "Carcinoma_Fractions", "Ecotypes", "discovery")

  all_H = NULL
  for(cell_type in key[,1]){
    #print(cell_type)
    # n_clusters = key[key[,1] == cell_type,2]
    classes = read.delim(file.path(states_dir, cell_type, 'state_abundances.txt'))
    rownames(classes) = paste0(cell_type, "_", rownames(classes))
    all_H = rbind(all_H, classes)
  }

  all_H = all_H[rownames(all_H) %in% Index$CellState,]

  all_H = t(all_H)
  colnames(all_H) = make.names(colnames(all_H))

  favor = Index$CellState[Index$Response == "Favor"]
  unfavor = Index$CellState[Index$Response == "Unfavor"]

  favor = apply(all_H[,favor],1,sum)
  unfavor = apply(all_H[,unfavor],1,sum)
  FavTME =  favor - unfavor
  FavTME = data.frame(FavTME = FavTME, row.names = rownames(all_H))

  colnames(all_H) = Index$Names[match(colnames(all_H),Index$CellState)]
  all_H = as.data.frame(all_H)

  return (list(FavTME = FavTME,
       CellStates = all_H))

}
