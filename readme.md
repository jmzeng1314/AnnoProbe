### 表达芯片数据分析伴侣

上周我们发布的四个R包基本上能进解决**表达芯片数据挖掘的88%的问题**，如下:

- [第三个万能芯片探针ID注释平台R包](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247492052&idx=1&sn=e640e0d9468f60616bddf9494117d0b8&chksm=9b4ba16fac3c287934abe0c6515f12af51d018ea2fa9436b2b9b0a1975c5a05318f2178448cd&scene=21#wechat_redirect)
- [第二个万能芯片探针ID注释平台R包](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247492037&idx=1&sn=f6a3a38cac4c20b5428f354803444ec4&chksm=9b4ba17eac3c286897898ac8c3cd42418600cf9b206ca3f9c97b04bcf8c49553fcd52ba891ea&scene=21#wechat_redirect)
- [第一个万能芯片探针ID注释平台R包](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247492027&idx=1&sn=9ae651dda053e3e3778d4557b141b6bd&chksm=9b4ba100ac3c28160ac6eaed032782c35218330420cf197d26fbff83602f35a556d25433fe0e&scene=21#wechat_redirect)
- [GEO数据库中国区镜像横空出世](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247492027&idx=2&sn=d8d1f2009aa9e50506c2259021b65e1a&chksm=9b4ba100ac3c281677f74a72910fa25662e9e7a888542e29867400458558593d3986a5f6e19c&scene=21#wechat_redirect)

有趣的是，因为这些包存储在GitHub，而且每个包自带的**数据是40~50M**，对很多在中国大陆的朋友来说， 几乎是不可能完成，所以我把这4个包整合成为了一个GitHub包（**AnnoProbe**）！总共不到5M，相信大家使用起来应该是很方便啦！

### 首先看GEO数据库下载镜像

```r
rm(list = ls())
library(AnnoProbe) 
suppressPackageStartupMessages(library(GEOquery)) 
gset=AnnoProbe::geoChina('GSE1009')
gset
```

需要理解一下下载的ExpressionSet对象，主要是表达矩阵和临床信息啦！

```r
# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
probes_expr=log2(probes_expr+1)
boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
```

### 然后对表达芯片的探针进行基因注释

这一个步骤也是在线下载我们的芯片注释信息

```r
## check GPL and annotate the probes to genes.
(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl)
head(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene )
head(genes_expr)
```

### 走limma的经典2组差异分析

```r
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
```

### 对差异分析结果进行一些检验

```R
## visualization
need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
deg_volcano(need_deg,1)
deg_volcano(need_deg,2)

deg_heatmap(DEG,genes_expr,group_list)
deg_heatmap(DEG,genes_expr,group_list,30)

check_diff_genes('PLCE1',genes_expr,group_list)
check_diff_genes('MPP6',genes_expr,group_list)
```

### 如果你做了GO/KEGG注释后也可以挑选基因集进行可视化

```r
# 假设我这里对hsa03410感兴趣
library(KEGGREST)
cg <- KEGGREST::keggGet("hsa03410")[[1]]$GENE
cg=as.character(sapply(cg[seq(2,length(cg),by=2)], function(x) strsplit(x,';')[[1]][1]))
check_diff_genes( cg ,genes_expr,group_list)
```

### 上面的代码你可以套用到任何一个表达芯片数据集

当然了，你需要有一点R语言基础知识啦，不然，上面的代码你不知道应该修改哪里。

马上试试看下面的数据集吧，我觉得蛮有意义的。

```
GSE1462
GSE18732
GSE20950
GSE21785
GSE26526
GSE32575
GSE43837
GSE474
GSE58979
GSE60291
GSE62832
GSE70529
GSE72158
```

