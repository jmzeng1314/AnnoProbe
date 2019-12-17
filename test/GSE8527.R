rm(list = ls())
library(AnnoProbe)
library(ggpubr)
suppressPackageStartupMessages(library(GEOquery))
getwd()
setwd('test/')
gset=AnnoProbe::geoChina('GSE8527')
gset
# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
probes_expr=log2(probes_expr+1)
boxplot(probes_expr,las=2)
library(limma)
probes_expr=normalizeBetweenArrays(probes_expr)
boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
group_list=ifelse(grepl('Control',phenoDat[,1]),'control','treat')

## check GPL and annotate the probes to genes.
(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl)
head(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene )
head(genes_expr)

# do DEG
## define the group
#group_list=factor(c(rep('Control',3),rep('Diabetes',3)))
group_list=ifelse(grepl('Control',phenoDat[,1]),'control','treat')

table(group_list)
library(limma)
design=model.matrix(~factor(group_list))
design
fit=lmFit(genes_expr,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
head(DEG)

## visualization
need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
deg_volcano(need_deg,1)
deg_volcano(need_deg,2)

deg_heatmap(DEG,genes_expr,group_list)
deg_heatmap(DEG,genes_expr,group_list,30)

check_diff_genes('NFKBIZ',genes_expr,group_list)
check_diff_genes('CCL2',genes_expr,group_list)
library(KEGGREST)
cg <- KEGGREST::keggGet("hsa03410")[[1]]$GENE
cg=as.character(sapply(cg[seq(2,length(cg),by=2)], function(x) strsplit(x,';')[[1]][1]))
check_diff_genes( cg ,genes_expr,group_list)






