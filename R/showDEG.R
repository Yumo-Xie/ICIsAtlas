#' @importFrom magrittr           %>%
#' @importFrom dplyr              filter select arrange slice_head mutate matches
#' @importFrom rlang              sym
#' @import                        ggplot2
#' @importFrom patchwork          wrap_plots
NULL

#' Plot meta-analysis differential expression for a gene between responders (Red) and non-responders (Blue).
#'
#' Visualizes log-fold changes and confidence intervals for a specified gene
#' across studies or within a chosen tumor type.
#'
#' @param gene character; gene symbol or ID to visualize.
#' @param Type character; one of \code{"Pan-cancer"}, \code{"Melanoma"}, \code{"Renal"}, or \code{"Urothelial"}.
#' @return Invisibly returns a \pkg{ggplot2} object; also prints the plot.
#' @examples
#' \dontrun{
#' showDEG(
#'   gene   = "PDCD1",
#'   Type   = "Pan-cancer"
#' )
#' }
#' @export
showDEG = function (gene, Type = "Pan-cancer")
{
  genecol = "ID"
  foldchangecol="logFC"
  llcol="CI.L"
  rlcol="CI.R"
  if (!exists("DEG_meta", inherits = TRUE)) {
    download_DEG_meta()
  }
  if(Type == "Pan-cancer"){
    meta = DEG_meta[[1]]
  }else if(Type == "Melanoma"){
    meta = DEG_meta[[2]]
  }else if(Type == "Renal"){
    meta = DEG_meta[[3]]
  }else if(Type == "Urothelial"){
    meta = DEG_meta[[4]]
  }else{
    stop(paste("Seems that", Type, "is not in the",
                            "tumor list"))
  }


  rem <- sremres <- merge(meta@metaresult, meta@input,
                          by = genecol) %>% dplyr::filter(!!rlang::sym(genecol) ==
                                                            gene)

  stds <- unique(unlist(regmatches(colnames(sremres), regexec("_\\d+$",
                                                              colnames(sremres)))))

  stds <- setNames(stds, meta@inputnames)
  edat <- Reduce(rbind, lapply(names(stds), function(sn) {
    std <- dplyr::select(sremres, dplyr::matches(paste0(genecol,
                                                        "|", stds[sn], "$")))
    colnames(std) <- gsub("_\\d+$", "", colnames(std))
    std[["group"]] <- sn
    std
  }))

  edat <- dplyr::select(edat, c(!!rlang::sym(genecol), !!rlang::sym(foldchangecol),
                                !!rlang::sym(llcol), !!rlang::sym(rlcol), group))
  edat <- edat %>%
    dplyr::mutate(
      Signconsistancy = NA_real_,Pval = NA_real_)
  sdat <- data.frame(genecol = unique(edat[[genecol]]), foldchangecol = sremres[["randomSummary"]],
                     llcol = sremres[["randomCi.lb"]], rlcol = sremres[["randomCi.ub"]],
                     group = "FoldChange summary",Signconsistancy = sremres[["signcon"]], Pval = sremres[["randomP"]])
  colnames(sdat) <- c(genecol, foldchangecol,  llcol, rlcol,
                      "group","Signconsistancy","Pval")
  dat <- rbind(edat, sdat)
  dat[["class"]] <- ifelse(grepl("summary", dat[["group"]]),
                           "FoldChange summary", "Study")
  sumfc <- dplyr::filter(dat, grepl("summary", class))[[foldchangecol]]
  maxfc <- max(dat[[rlcol]])
  minfc <- min(dat[[llcol]])
  if (sumfc > 0) {
    sumcol <- "#E41A1C"
    minlim <- -maxfc
    maxlim <- maxfc
  }
  else {
    sumcol <- "#377EB8"
    minlim <- minfc
    maxlim <- -minfc
  }
  dat[!is.na(dat$Pval),"Pval"] = paste0("P-value =", " ", round(dat[!is.na(dat$Pval),"Pval"],digits = 3))
  dat[!is.na(dat$Signconsistancy),"Signconsistancy"] = paste0("Consistancy =", " ",dat[!is.na(dat$Signconsistancy),"Signconsistancy"])
  gg <- ggplot(dat, aes(x = group, y = !!rlang::sym(foldchangecol),
                        color = class)) + geom_point() + geom_errorbar(aes(ymin = !!rlang::sym(llcol),
                                                                           ymax = !!rlang::sym(rlcol), width = 0.1, color = class)) +
    geom_text(data = subset(dat, !is.na(Signconsistancy) & !is.na(Pval)), aes(x = group,y = maxlim*0.8,
                    label = paste0("logFC =", " ",round(logFC,2))),color = sumcol,size = 5,inherit.aes = FALSE, vjust = -3) +
    geom_text(data = subset(dat, !is.na(Signconsistancy) & !is.na(Pval)), aes(x = group,y = maxlim*0.8,
                    label = Pval),color = sumcol,size = 5,inherit.aes = FALSE, vjust = -1.25) +
    geom_text(data = subset(dat, !is.na(Signconsistancy) & !is.na(Pval)), aes(x = group,y = maxlim*0.8,
                    label = Signconsistancy),color = sumcol,size = 5,inherit.aes = FALSE) +
    scale_color_manual(values = c(sumcol, "#bdbdbd")) + scale_x_discrete(limits = rev(dat[["group"]])) +
    theme_classic() + ggtitle(unique(edat[[genecol]])) + theme(axis.title.y = element_blank()) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.05,
               color = "#969696") + geom_hline(yintercept = sumfc,
                                               linetype = "dashed", size = 0.1, color = sumcol) + theme(legend.position = "none") +
    scale_y_continuous(limits = c(minlim, maxlim)) + coord_flip()


  if (!exists("DEG_meta", envir = .GlobalEnv)) {
    return(DEG_meta)
  }
   return(gg)
}

#' Plot pathway enrichment results
#'
#' Generates enrichment results for specified pathways between responders (Red) and non-responders (Blue), or automatically
#' selects the top \code{Top} pathways if \code{Pathway = NULL}.
#'
#' @param Pathway character vector of pathway names; if \code{NULL}, uses \code{Top}.
#' @param Top integer; number of top (positive) or bottom (negative) pathways to select when \code{Pathway = NULL}.
#' @param Type character; tumor type, same options as in \code{showDEG}.
#' @param stack logical; if \code{TRUE}, returns a combined patchwork object; otherwise prints individual plots.
#' @return If \code{stack = TRUE}, returns a \pkg{patchwork} object; else invisibly \code{NULL}.
#' @examples
#' \dontrun{
#' showDEP(
#'   Pathway = "GOBP_NATURAL_KILLER_CELL_CHEMOTAXIS",
#'   Type    = "Pan-cancer",
#'   stack   = F
#' )
#' }
#' @export
showDEP = function (Pathway = NULL, Top = 5, Type = "Pan-cancer",stack = F)
{
  genecol = "var"
  foldchangecol="logFC"
  llcol="CI.L"
  rlcol="CI.R"
  if(Type == "Pan-cancer"){
    meta = GO_meta[[1]]
  }else if(Type == "Melanoma"){
    meta = GO_meta[[2]]
  }else if(Type == "Renal"){
    meta = GO_meta[[3]]
  }else if(Type == "Urothelial"){
    meta = GO_meta[[4]]
  }else{
    stop(paste("Seems that", Type, "is not in the",
               "tumor list"))
  }

  if(is.null(Pathway) & Top > 0){
    Pathway = meta@metaresult %>% dplyr::arrange(desc(randomSummary)) %>% dplyr::slice_head(n = abs(Top)) %>% dplyr::pull(var)
  }else if(is.null(Pathway) & Top < 0){
    Pathway = meta@metaresult %>% dplyr::arrange(randomSummary) %>% dplyr::slice_head(n = abs(Top)) %>% dplyr::pull(var)
  }

  plot_list <- vector("list", length = length(Pathway))
  for (i in 1:length(Pathway)) {

  rem <- sremres <- merge(meta@metaresult, meta@input,
                          by = genecol) %>% dplyr::filter(!!rlang::sym(genecol) ==
                                                            Pathway[i])

  stds <- unique(unlist(regmatches(colnames(sremres), regexec("_\\d+$",
                                                              colnames(sremres)))))
  stds <- setNames(stds, meta@inputnames)
  edat <- Reduce(rbind, lapply(names(stds), function(sn) {
    std <- dplyr::select(sremres, dplyr::matches(paste0(genecol,
                                                        "|", stds[sn], "$")))
    colnames(std) <- gsub("_\\d+$", "", colnames(std))
    std[["group"]] <- sn
    std
  }))

  edat <- dplyr::select(edat, c(!!rlang::sym(genecol), !!rlang::sym(foldchangecol),
                                !!rlang::sym(llcol), !!rlang::sym(rlcol), group))
  edat <- edat %>%
    dplyr::mutate(
      Signconsistancy = NA_real_,Pval = NA_real_)
  sdat <- data.frame(genecol = unique(edat[[genecol]]), foldchangecol = sremres[["randomSummary"]],
                     llcol = sremres[["randomCi.lb"]], rlcol = sremres[["randomCi.ub"]],
                     group = "FoldChange summary",Signconsistancy = sremres[["signcon"]], Pval = sremres[["randomP"]])
  colnames(sdat) <- c(genecol, foldchangecol,  llcol, rlcol,
                      "group","Signconsistancy","Pval")
  dat <- rbind(edat, sdat)
  dat[["class"]] <- ifelse(grepl("summary", dat[["group"]]),
                           "FoldChange summary", "Study")
  sumfc <- dplyr::filter(dat, grepl("summary", class))[[foldchangecol]]
  maxfc <- max(dat[[rlcol]])
  minfc <- min(dat[[llcol]])
  if (sumfc > 0) {
    sumcol <- "#E41A1C"
    minlim <- -maxfc
    maxlim <- maxfc
  }
  else {
    sumcol <- "#377EB8"
    minlim <- minfc
    maxlim <- -minfc
  }
  dat[!is.na(dat$Pval),"Pval"] = paste0("P-value =", " ", round(dat[!is.na(dat$Pval),"Pval"],digits = 3))
  dat[!is.na(dat$Signconsistancy),"Signconsistancy"] = paste0("Consistancy =", " ",dat[!is.na(dat$Signconsistancy),"Signconsistancy"])
  plot_list[[i]] <- ggplot(dat, aes(x = group, y = !!rlang::sym(foldchangecol),
                        color = class)) + geom_point() + geom_errorbar(aes(ymin = !!rlang::sym(llcol),
                                                                           ymax = !!rlang::sym(rlcol), width = 0.1, color = class)) +
    geom_text(data = subset(dat, !is.na(Signconsistancy) & !is.na(Pval)), aes(x = group,y = maxlim*0.7,
                                                                              label = paste0("logFC =", " ",round(logFC,2))),color = sumcol,size = 5,inherit.aes = FALSE, vjust = -3) +
    geom_text(data = subset(dat, !is.na(Signconsistancy) & !is.na(Pval)), aes(x = group,y = maxlim*0.7,
                                                                              label = Pval),color = sumcol,size = 5,inherit.aes = FALSE, vjust = -1.25) +
    geom_text(data = subset(dat, !is.na(Signconsistancy) & !is.na(Pval)), aes(x = group,y = maxlim*0.7,
                                                                              label = Signconsistancy),color = sumcol,size = 5,inherit.aes = FALSE) +
    scale_color_manual(values = c(sumcol, "#bdbdbd")) + scale_x_discrete(limits = rev(dat[["group"]])) +
    theme_classic() + ggtitle(unique(edat[[genecol]])) + theme(axis.title.y = element_blank()) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.05,
               color = "#969696") + geom_hline(yintercept = sumfc,
                                               linetype = "dashed", size = 0.1, color = sumcol) + theme(legend.position = "none") +
    scale_y_continuous(limits = c(minlim, maxlim)) + coord_flip()

   if(stack == F){
     print(plot_list[[i]])
   }
  }
  if(stack == T){
    gg = patchwork::wrap_plots(plot_list, nrow = i)
    return(gg)
  }

}

