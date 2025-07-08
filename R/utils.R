#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomicRanges   reduce width
#' @import        org.Hs.eg.db
#' @import        TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom edgeR           DGEList calcNormFactors rpkm
NULL

#' @keywords internal
#' @noRd
cal_rpkm = function(expr){
  ex <- exonsBy(tbdx, by="gene")
  gene_length <- sum(width(reduce(ex)))
  gene_length <- as.numeric(gene_length)
  names(gene_length) <- names(ex)
  symbols <- rownames(expr)
  map_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("ENTREZID")
  )
  map_df <- map_df[!is.na(map_df$ENTREZID), ]
  map_df$Length <- gene_length[map_df$ENTREZID]
  map_df <- map_df[ !duplicated(map_df$SYMBOL), ]
  names(gene.length) <- map_df$SYMBOL
  gene.length <- gene.length[symbols]

  dge <- DGEList(counts=expr, genes = data.frame(Symbol=symbols))
  dge$genes$Length <- gene.length
  dge = calcNormFactors(dge)
  expr <- rpkm(dge, gene.length="Length")

  return(expr)

}

#' @keywords internal
#' @noRd
mypval=function(x){

  ifelse(x < 0.01,
         ifelse(x < 0.001,"p < 0.001",paste("p =",round(x,digits = 3),sep = " ")),
         paste("p =",round(x,digits = 2),sep = " "))
}
