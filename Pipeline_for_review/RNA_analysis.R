library(limma)
library(DESeq2)
library(edgeR)
library(MetaVolcanoR)
library(dplyr)
library(HGNChelper)

update_gene_symbols = function(mat, type = c("sum", "average")) {
  gene_list = rownames(mat)
  corrected = checkGeneSymbols(gene_list, unmapped.as.na = TRUE)
  new_symbols = ifelse(is.na(corrected$Suggested.Symbol), 
                        corrected$x, 
                        corrected$Suggested.Symbol)
  
  if (type == "sum") {
    mat = rowsum(mat, group = new_symbols)
  } else {
    mat = limma::avereps(mat, ID = new_symbols)
  }
  return(mat)
}


# Count data =============

DEGlist=list()
for (i in 1:length(count_list)) {
  mat=count_list[[i]]$mat
  mat=update_gene_symbols(mat)
  if (max(mat, na.rm = TRUE) > 50) {
    mat = log2(mat + 1)
  }
  
  condition = Data_list[[i]]$res
  expfactor = factor(condition, levels = c('nonresponder', 'responder'))
  
  des = model.matrix(~0 + expfactor)
  colnames(des) = c('nonresponder', 'responder')
  
  dgelist = DGEList(counts = mat,group = expfactor)
  keep = rowSums(cpm(dgelist) > 1 ) >= 2
  dgelist = dgelist[keep, , keep.lib.sizes = FALSE]
  
  
  dgelist_norm = calcNormFactors(dgelist, method = 'TMM')
  
  logCPM = voom(dgelist_norm, des, plot = TRUE)
  
  contrast.matrix=makeContrasts("responder-nonresponder",
                                 levels = des)
  
  fit = lmFit(logCPM, des)
  fit2 = contrasts.fit(fit, contrast.matrix) 
  fit2 = eBayes(fit2)
  diffexp = topTable(fit2, coef="responder-nonresponder", n=Inf, confint=0.95) %>% na.omit()
  diffexp$ID = rownames(diffexp)
  DEGlist[[i]]=diffexp
  names(DEGlist)[i]=names(count_list[i])
}


# TPM and array--------------------------------

library(genefilter)

DEGlist1=list()
for (i in 1:length(Data_list)) {
  mat=Data_list[[i]]$mat
  mat=update_gene_symbols(mat)
  
  expfactor = Data_list[[i]]$res
  expfactor = factor(expfactor, levels = c('nonresponder', 'responder'))
  des = model.matrix(~0 + expfactor)
  colnames(des) = c('nonresponder', 'responder')
  
  contrast.matrix=makeContrasts("responder-nonresponder",
                                 levels = des)
  print(dim(mat))
  
  if(max(mat)>50){
    mat=log2(mat+1)
    message("Perform log transform for ", names(Data_list)[i])
  }
  
  mat1 = varFilter(mat,var.cutoff=0.3)
  print(dim(mat1))
  
  fit = lmFit(mat1, des)
  fit2 = contrasts.fit(fit, contrast.matrix) 
  fit2 = eBayes(fit2, trend = TRUE)
  print(plotSA(fit2, main = "Mean-variance trend") )
  diffexp1 = topTable(fit2, coef="responder-nonresponder", n=Inf, confint=0.95) %>% na.omit()
  diffexp1$ID = rownames(diffexp1)
  DEGlist1[[i]]=diffexp1
  names(DEGlist1)[i]=names(Data_list[i])
}



# Meta =========

metalist= c(DEGlist,DEGlist1)

meta_degs_rem = rem_mv(diffexp=metalist,
                        pcriteria="P.Value",
                        foldchangecol='logFC', 
                        genenamecol='ID',
                        geneidcol=NULL,
                        collaps=T,
                        llcol='CI.L',
                        rlcol='CI.R',
                        vcol=NULL, 
                        cvar=TRUE,
                        metathr=0.01,
                        jobname="MetaVolcano",
                        outputfolder=".", 
                        draw='PDF',
                        ncores=48)

result=meta_degs_rem@metaresult
result=result[result$randomP < 0.05,]
result=result[result$signcon > 4 | result$signcon < -4,] 
result=result[order(result$randomSummary,decreasing = T),]

