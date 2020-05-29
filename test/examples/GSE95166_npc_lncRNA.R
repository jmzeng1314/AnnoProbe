rm(list = ls())
library(AnnoProbe)
library(ggpubr)
suppressPackageStartupMessages(library(GEOquery))
getwd()
setwd('test/')
## download GSE95166 data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95166
eSet=getGEO('GSE95166', destdir=".", AnnotGPL = F, getGPL = F)[[1]]
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
# https://www.ncbi.nlm.nih.gov/pubmed/31430288

group_list=factor(c(rep('npc',4),rep('normal',4)))
table(group_list)
eSet@annotation
# GPL15314	Arraystar Human LncRNA microarray V2.0 (Agilent_033010 Probe Name version)

GPL=eSet@annotation
probes_anno <- idmap(GPL,type = 'pipe')
head(probes_anno)
genes_expr <- filterEM(probes_expr,probe2gene )
head(genes_expr)

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
# We observed that 2107 lncRNAs were upregulated
# while 2090 lncRNAs were downregulated by more than 2-fold,
# NKILA among these downregulated lncRNAs (Fig 1A, GSE95166).


## visualization
need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
deg_volcano(need_deg,1)
deg_volcano(need_deg,2)

deg_heatmap(DEG,genes_expr,group_list)
deg_heatmap(DEG,genes_expr,group_list,30)

check_diff_genes('NKILA',genes_expr,group_list)
cg=c('CR612603','NKILA',
     'PR11-392A14.3',
     'AC013474.3',
     'RP3-351K20.3',
     'RP11-45P22.1',
     'CR624187',
     'KRT16P1',
     'RP11-215C7.2',
     'AL359062',
     'AC004797.1',
     'AK123324',
     'ZNF295-AS1',  # 'C21orf121',
     'SMIM6', # 'AK131023',
     'BC022056',
     'FAM47E',#'LincRNA-FAM47E',
     'RP11-392A23.4',
     'BC041457',
     'AK095282',
     'GCNT3',#'lincrns-GCNT3',
     'ABCC6P2')
check_diff_genes( cg ,genes_expr,group_list)

tail(sort(table(annoGene(need_deg$symbols,'SYMBOL')[,2])))




