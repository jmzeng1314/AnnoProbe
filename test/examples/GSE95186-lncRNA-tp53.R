rm(list = ls())
library(AnnoProbe)
library(ggpubr)
suppressPackageStartupMessages(library(GEOquery))
getwd()
setwd('test/')
## download GSE95186 data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95186
# eSet=getGEO('GSE95186', destdir=".", AnnotGPL = F, getGPL = F)[[1]]
# 使用 GEO 中国区镜像进行加速

gset=AnnoProbe::geoChina('GSE95186')
gset
# check the ExpressionSet
eSet=gset[[1]]
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
# https://www.ncbi.nlm.nih.gov/pubmed/31430288

group_list=factor(c(rep('treat',3),rep('untreat',3)))
table(group_list)
eSet@annotation
# GPL15314	Arraystar Human LncRNA microarray V2.0 (Agilent_033010 Probe Name version)

GPL=eSet@annotation
# 选择 pipe 获取的是 冗余注释，也就是说一个探针很有可能会对应多个基因。
probes_anno <- idmap(GPL,type = 'pipe')
head(probes_anno)
probes_anno=probes_anno[probes_anno$probe_id %in% rownames(probes_expr),]
# 只需要表达矩阵里面有的探针的注释即可

length(unique(probes_anno$probe_id))
anno=annoGene(probes_anno$symbol,'SYMBOL')
head(anno)
pcs=probes_anno[probes_anno$symbol %in% anno[anno$biotypes=='protein_coding',1],]
nons=probes_anno[probes_anno$symbol %in% anno[anno$biotypes !='protein_coding',1],]
# 可以首先把探针拆分成为 protein_coding 与否
length(unique(pcs$probe_id))
length(unique(nons$probe_id))
# 可以看到仍然是有探针会被注释到多个基因，这个时候

pcs_expr <- probes_expr[pcs$probe_id,]
nons_expr <- probes_expr[nons$probe_id,]
boxplot(pcs_expr,las=2)
boxplot(nons_expr,las=2)
# 很容易看出来，非编码的这些基因的平均表达量，是低于编码的。

## 首先对编码基因的表达矩阵做差异分析
genes_expr <- filterEM(pcs_expr,pcs )
if(T){
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
        ## visualization
        need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
        deg_volcano(need_deg,1)
        deg_volcano(need_deg,2)

        deg_heatmap(DEG,genes_expr,group_list)
        deg_heatmap(DEG,genes_expr,group_list,30)
}

### 然后对非编码的基因的表达矩阵做差异分析
genes_expr <- filterEM(nons_expr,nons )
if(T){
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
        ## visualization
        need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
        deg_volcano(need_deg,1)
        deg_volcano(need_deg,2)

        deg_heatmap(DEG,genes_expr,group_list)
        deg_heatmap(DEG,genes_expr,group_list,30)
}





