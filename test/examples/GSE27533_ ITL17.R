rm(list = ls())
library(AnnoProbe)
library(ggpubr)
suppressPackageStartupMessages(library(GEOquery))
getwd()
setwd('test/')
## download GSE27533 data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27533
eSet=getGEO('GSE27533', destdir=".", AnnotGPL = F, getGPL = F)[[1]]
# GSE27533=eSet;usethis::use_data(GSE27533)
## expgenes_exprsion matrix
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
# https://www.ncbi.nlm.nih.gov/pubmed/31430288

group_list=factor(c(rep('Control',3),rep('IL17A',3)))
table(group_list)

genes_expr=probes_expr

library("FactoMineR")
library("factoextra")
dat.pca <- PCA(t(genes_expr) , graph = FALSE)
dat.pca
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             legend.title = "Groups"
)
library(limma)
design=model.matrix(~factor(group_list))
design
fit=lmFit(genes_expr,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
head(DEG)

eSet@annotation
# 	Illumina HumanHT-12 V4.0 expression beadchip
GPL=eSet@annotation
probes_anno <- idmap(GPL,type = 'pipe')
head(probes_anno)
DEG$probe_id=rownames(DEG)
DEG=merge(DEG,probes_anno,by='probe_id')

## visualization
need_deg=data.frame(symbols=DEG$symbol, logFC=DEG$logFC, p=DEG$P.Value)
deg_volcano(need_deg,1)
deg_volcano(need_deg,2)

deg_heatmap(DEG,genes_expr,group_list)
deg_heatmap(DEG,genes_expr,group_list,30)

# ILMN_1674038 CTSD
probes_expr['ILMN_1674038',]
genes_expr['CTSD',]
check_diff_genes('GAPDH',genes_expr,group_list)
cg=c('CTSD','CAPNS1','KLK5','CAPNS1','CTSL1','CTSB')
check_diff_genes( cg ,genes_expr,group_list)

tail(sort(table(annoGene(need_deg$symbols,'SYMBOL')[,2])))

