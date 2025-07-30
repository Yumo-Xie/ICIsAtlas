#' @importFrom magrittr           %>%
#' @importFrom dplyr              filter mutate if_else
#' @importFrom rlang              sym
#' @importFrom survminer          surv_cutpoint surv_categorize ggsurvplot surv_pvalue
#' @importFrom survival           Surv survfit
#' @importFrom crayon             green red
NULL

#' Plot Kaplan–Meier survival curves
#'
#' KM curves (and risk tables) stratified by a chosen metric,
#' either per-study or pooled across all samples.
#'
#' @param survival character; outcome, either \code{"OS"} or \code{"PFS"}.
#' @param type character; stratification variable (e.g.\ \code{"FavTME"} or cell-state name).
#' @param cut character; cutpoint method, \code{"median"} or \code{"best"}.
#' @param by.study logical; if \code{TRUE}, plots curves separately per study; else pooled.
#' @param available logical; if \code{TRUE}, prints available input types and exits invisibly.
#' @return Returns plots.
#' @examples
#' \dontrun{
#' # Using available = T to see the available input types -----------------
#'
#' survplot(available = T)
#'
#' # Plot survival curves ----------------------
#' survplot(
#'   survival = "PFS",
#'   type     = "FavTME",
#'   cut      = "median",
#'   by.study = T
#' )
#' }
#' @export
survplot = function(survival = c("OS","PFS"), type = "FavTME", cut = c("median","best"), by.study = T, available = F)
  {
  survival = match.arg(survival)
  cut = match.arg(cut)
  if(available == T){
    cat(crayon::green$italic("Available input types include..."), "\n")
    cat(crayon::red("FavTME"), "\n")
    cat(crayon::green$italic("Cellstates in EcoTyper"), "\n")
    print(Index,n = 69)
    cat(crayon::green$italic("Cells in LM22"), "\n")
    print(LM22)
    cat(crayon::green$italic("To perform survival analysis, set available to FALSE"), "\n")
    return(invisible(NULL))
  }

  if(!(type %in% Index$Names) & !(type %in% c("FavTME")) & !(type %in% Index$CellState) & !(type %in% LM22)){
    stop(paste("Seems that", type, "is not in the",
               "cell/FavTME lists"))
  }

  if(type %in% Index$CellState){
    type = Index$Names[match(type,Index$CellState)]
  }else{
    type_input = type
    type = make.names(type)
  }

  if(!(cut %in% c("median","best"))){
    stop("Incorrect cutpoint")
  }
  if(survival %in% c("OS","PFS")){
    survival_status = paste0(survival,"_status")
  }else{
    stop("Incorrect survival input")
  }

  if(by.study == T){
    # if(survival == "OS"){
    #   data = OS_list
    # }else{
    #   data = PFS_list
    # }

    data = Survival_list

    for (i in 1:length(data)) {
      rt=data[[i]]

      ylab = ifelse(survival == "PFS", "Progression-free Survial (%)","Overall Survial (%)")
      ptitle = paste0(unique(rt$Study),"   ",type_input)

      if(all(is.na(rt[,survival]))){next} else if(cut == "median"){
          rt = rt %>% dplyr::mutate(!!rlang::sym(type):= if_else(
            !!rlang::sym(type) > median(!!rlang::sym(type)),
            "high","low"))
        }else{
          rt = survminer::surv_cutpoint(rt,
                                  time = survival,
                                  event = survival_status,
                                  variables = type)
          rt = survminer::surv_categorize(rt)
        }

        fit <- survfit(as.formula(paste0(
          "Surv(", survival, ", ", survival_status, ") ~ ", type)), data = rt)

        fit$call$formula = eval(fit$call$formula)

        if(max(rt[,survival],na.rm = T) < 10){breaks=3}else{if(max(rt[,survival],na.rm = T) > 50){breaks=12}else{breaks=6}}

        b=ggsurvplot(fit,
                     data = rt,
                     conf.int = FALSE,
                     fun = "pct",
                     # pval = TRUE,
                     censor.shape = "|",
                     censor.size=2,
                     risk.table = TRUE,
                     xlim=c(0,ceiling(max(rt[,survival],na.rm = T)/breaks)*breaks),
                     # ylim=c(50,100),
                     break.x.by=breaks,
                     size = 1,
                     palette = c("#BC3C29FF","#0072B5FF"),
                     title= ptitle,
                     ylab= ylab,
                     xlab="Time (months)",
                     legend.labs = c("High", "Low"),
                     legend=c(0.8,1),
                     legend.title= "",
                     pval=mypval(surv_pvalue(fit,data = rt)$pval),
                     pval.size=4,
                     pval.coord=c(ceiling(max(rt[,survival],na.rm = T)/breaks)*breaks*0.8,85),
                     font.tickslab=16,
                     font.y=16,
                     font.x=16,
                     fontsize=6,
                     tables.height=0.3,
                     font.table.x=6,
        )
        b$table$theme$axis.title.x$size=16
        b$table$theme$axis.text.x$size=16
        b$table$theme$axis.text.y$size=16
        print(b)

      }

  }else if(by.study == F){


    rt = Survival_data
    rt = rt %>% filter(!is.na(!!rlang::sym(survival)))

    if(cut == "median"){
      rt = rt %>% dplyr::mutate(!!rlang::sym(type):= if_else(
        !!rlang::sym(type) > median(!!rlang::sym(type)),
                                    "high",
                                    "low"
        ))
    }else{
      rt = survminer::surv_cutpoint(rt,
                                    time = survival,
                                    event = survival_status,
                                    variables = type)
      rt = survminer::surv_categorize(rt)
    }
    fml_text <- paste0(
      "Surv(", survival, ", ", survival_status, ") ~ ", type)

    fml <- as.formula(fml_text)
    fit <- survfit(fml, data = rt)
    fit$call$formula = eval(fit$call$formula)

    if(max(rt[,survival],na.rm = T) < 10){breaks=3}else{if(max(rt[,survival],na.rm = T) > 50){breaks=12}else{breaks=6}}

    ylab = ifelse(survival == "PFS", "Progression-free Survial (%)","Overall Survial (%)")

    b=ggsurvplot(fit,
                 data = rt,
                 conf.int = FALSE,
                 fun = "pct",
                 # pval = TRUE,
                 censor.shape = "|",
                 censor.size=2,
                 risk.table = TRUE,
                 xlim=c(0,ceiling(max(rt[,survival],na.rm = T)/breaks)*breaks),
                 # ylim=c(50,100),
                 break.x.by=breaks,
                 size = 1,
                 palette = c("#BC3C29FF","#0072B5FF"),
                 title= type_input,
                 ylab= ylab,
                 xlab="Time (months)",
                 legend.labs = c("High", "Low"),
                 legend=c(0.8,1),
                 legend.title= "",
                 pval=mypval(surv_pvalue(fit,data = rt)$pval),
                 pval.size=4,
                 pval.coord=c(ceiling(max(rt[,survival],na.rm = T)/breaks)*breaks*0.8,85),
                 font.tickslab=16,
                 font.y=16,
                 font.x=16,
                 fontsize=6,
                 tables.height=0.3,
                 font.table.x=6,
    )
    b$table$theme$axis.title.x$size=16
    b$table$theme$axis.text.x$size=16
    b$table$theme$axis.text.y$size=16
    print(b)

  }

}

