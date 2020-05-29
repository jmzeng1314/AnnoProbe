# Agilent-028004 SurePrint G3 Human GE 8x60K Microarray (Feature Number version)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13607
rm(list = ls())
library(AnnoProbe)

ids2=idmap(gpl='GPL13607',type = 'soft')
load('GPL13607.Rdata')
library(GEOquery)
colnames(Table(gpl))
head(Table(gpl))
tmp=Table(gpl)
write.table(tmp,file = 'GPL13607.txt',sep = '\t',quote=F,row.names = F)
ids=Table(gpl)[,c(1,6)]
ids2=ids
head(ids2)

ids3=idmap(gpl='GPL13607',type = 'pipe')
head(ids3)
length(unique(ids2[,1]))
length(unique(ids3[,1]))

head(ids2)
## 首先查看，哪些ids2里面的注释，并不在ids3，因为ids3是一对多的映射关系。
table( paste(ids2[,1],ids2[,2],sep = '_') %in% paste(ids3[,1],ids3[,2],sep = '_') )
## 一致的注释有25054个探针，不一致的有16054，差异很显著，值得深思。
ids2_not_ids3=ids2[!paste(ids2[,1],ids2[,2],sep = '_') %in% paste(ids3[,1],ids3[,2],sep = '_'),]
table(nchar(ids2_not_ids3$symbol)>0)
# 可以看到，其中5882个探针，在ids2里面的是注释不到基因的。
# 另外的10172探针是注释 冲突的。
m2=merge(ids2_not_ids3,ids3,by.x="ID",by.y='probe_id')
table(nchar(m2$symbol.x)>0)
tmp1=annoGene(m2[,2],'SYMBOL')
tmp2=annoGene(m2[,3],'SYMBOL')
tail(sort(table(tmp2$biotypes)))




# 其中一个不能注释，这样的冲突，是因为基因名过时。
m2_not=m2[!m2[,2] %in% tmp1[,1],]
head(m2_not)
an=annoGene(m2_not$symbol,'SYMBOL','human')
m2_not=merge(m2_not,an,by.x = 'symbol',by.y='SYMBOL')
# 都能注释，但是确实是冲突
m2_ok=m2[m2[,2] %in% tmp1[,1],]
head(m2_ok)
an=annoGene(m2_ok$symbol,'SYMBOL','human')
m2_ok=merge(m2_ok,an,by.x = 'symbol',by.y='SYMBOL')









