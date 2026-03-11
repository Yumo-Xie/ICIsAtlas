library(MetaVolcanoR)
library(GSVA)
library(limma)
library(GSEABase)

gmtFile="c5.go.bp.v7.4.symbols.gmt" 
geneSet=getGmt(file.path(homedir,gmtFile), geneIdType=SymbolIdentifier())

metalist = list()
for (i in 1:length(lastlst)) {
  rt=lastlst[[i]]$matrix
  ssgseaScore=gsva(rt, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE, ssgsea.norm=TRUE, parallel.sz=6L)
  
  res=lastlst[[i]]$res
  
  design=model.matrix(~res)
  colnames(design)[2]="resVSnonres"
  
  fit <- lmFit(ssgseaScore, design)
  fit <- eBayes(fit)
  sig.total=topTable(fit, coef="resVSnonres",confint=T,number=7500)
  metalist[[i]]=sig.total
  names(metalist)[[i]]=names(lastlst[i])
}

metalist1=list()
for (i in 1:length(metalist)) {
  sig=metalist[[i]]
  sig$var=rownames(sig)
  metalist1[[i]]=sig
  names(metalist1)[[i]]=names(metalist[i])
}


meta_degs_rem2 = rem_mv(diffexp=metalist1,
                         pcriteria="adj.P.Val",
                         foldchangecol='logFC', 
                         genenamecol='var',
                         geneidcol=NULL,
                         collaps=T,
                         llcol='CI.L',
                         rlcol='CI.R',
                         vcol=NULL, 
                         cvar=TRUE,
                         metathr=0.01,
                         jobname="MetaVolcano",
                         outputfolder=".", 
                         draw='HTML',
                         ncores=6)

meta_degs_rem2@MetaVolcano
result2=meta_degs_rem2@metaresult

result2=result2[result2$randomP < 0.05,]
result2=result2[result2$signcon > 7 | result2$signcon < -7,]
