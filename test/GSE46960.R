rm(list = ls())
library(AnnoProbe)
library(ggpubr)
suppressPackageStartupMessages(library(GEOquery))
getwd()
setwd('test/')
gset=AnnoProbe::geoChina('GSE46960')
gset
# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr[,1:4],las=2)
# 表达矩阵被zscore，是作者的问题
# 可以从cel文件开始，重新处理
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])

## check GPL and annotate the probes to genes.
(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl)
# 因为默认的bioconductor提供的注释只有两万不到
# 所以换成soft的注释，接近3万个
probe2gene=idmap(gpl,type = 'soft')
head(probe2gene)
dim(probes_expr)
dim(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene )
head(genes_expr)
## 后面的分析就没有意义了，因为表达矩阵被zscore了。

# do DEG
## define the group
group_list=factor(c(rep('Control',3),rep('Diabetes',3)))
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

check_diff_genes('PLCE1',genes_expr,group_list)
check_diff_genes('MPP6',genes_expr,group_list)
library(KEGGREST)
cg <- KEGGREST::keggGet("hsa03410")[[1]]$GENE
cg=as.character(sapply(cg[seq(2,length(cg),by=2)], function(x) strsplit(x,';')[[1]][1]))
check_diff_genes( cg ,genes_expr,group_list)






