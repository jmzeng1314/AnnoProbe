rm(list = ls())
library(AnnoProbe)
ids1=idmap('GPL6947',type = 'bioc')
head(ids1)
ids2=idmap('GPL6947',type = 'soft')
head(ids2)
ids3=idmap('GPL6947',type = 'pipe')
head(ids3)
length(unique(ids1[,1]))
length(unique(ids2[,1]))
length(unique(ids3[,1]))


ids2_filter=ids2[nchar(ids2[,2])>1,]
m1=merge(ids1,ids2_filter,by.x = 'probe_id',by.y = 'ID')
table(m1[,2]==m1[,3])
m1_not=m1[m1[,2] !=m1[,3],]
tmp1=annoGene(m1_not[,2],'SYMBOL')
tmp2=annoGene(m1_not[,3],'SYMBOL')

table( paste(ids1[,1],ids1[,2],sep = '_') %in% paste(ids3[,1],ids3[,2],sep = '_') )
ids1_not_ids3=ids1[!paste(ids1[,1],ids1[,2],sep = '_') %in% paste(ids3[,1],ids3[,2],sep = '_'),]
m2=merge(ids1_not_ids3,ids3,by='probe_id')
tmp1=annoGene(m2[,2],'SYMBOL')
tmp2=annoGene(m2[,3],'SYMBOL')


m2_not=m2[!m2[,2] %in% tmp1[,1],]
m2_ok=m2[m2[,2] %in% tmp1[,1],]



rm(list = ls())
library(AnnoProbe)
ids1=idmap('GPL10558',type = 'bioc')
head(ids1)
ids2=idmap('GPL10558',type = 'soft')
head(ids2)
ids3=idmap('GPL10558',type = 'pipe')
head(ids3)
length(unique(ids1[,1]))
length(unique(ids2[,1]))
length(unique(ids3[,1]))


ids2_filter=ids2[nchar(ids2[,2])>1,]
m1=merge(ids1,ids2_filter,by.x = 'probe_id',by.y = 'ID')
table(m1[,2]==m1[,3])
m1_not=m1[m1[,2] !=m1[,3],]
tmp1=annoGene(m1_not[,2],'SYMBOL')
tmp2=annoGene(m1_not[,3],'SYMBOL')

table( paste(ids1[,1],ids1[,2],sep = '_') %in% paste(ids3[,1],ids3[,2],sep = '_') )
ids1_not_ids3=ids1[!paste(ids1[,1],ids1[,2],sep = '_') %in% paste(ids3[,1],ids3[,2],sep = '_'),]
m2=merge(ids1_not_ids3,ids3,by='probe_id')
tmp1=annoGene(m2[,2],'SYMBOL')
tmp2=annoGene(m2[,3],'SYMBOL')


m2_not=m2[!m2[,2] %in% tmp1[,1],]
m2_ok=m2[m2[,2] %in% tmp1[,1],]


head(ids2)
## 首先查看，哪些ids2里面的注释，并不在ids3，因为ids3是一对多的映射关系。
table( paste(ids2[,1],ids2[,2],sep = '_') %in% paste(ids3[,1],ids3[,2],sep = '_') )
## 一致的注释有 24738 个探针，不一致的有 23369 ，差异很显著，值得深思。
ids2_not_ids3=ids2[!paste(ids2[,1],ids2[,2],sep = '_') %in% paste(ids3[,1],ids3[,2],sep = '_'),]
table(nchar(ids2_not_ids3$symbol)>0)
# 可以看到，其中 3270 个探针，在ids2里面的是注释不到基因的。
# 另外的 20099 探针是注释 冲突的。
m2=merge(ids2_not_ids3,ids3,by.x="ID",by.y='probe_id')
table(nchar(m2$symbol.x)>0)
tmp1=annoGene(m2[,2],'SYMBOL')
tmp2=annoGene(m2[,3],'SYMBOL')
tail(sort(table(tmp2$biotypes)))
# 其中一种不能注释，这样的冲突，是因为基因名过时。
m2_not=m2[!m2[,2] %in% tmp1[,1],]
# 我们关心的是探针的冲突：16529个探针。
length(unique(m2_not$ID))
# 都能注释，但是确实是冲突
m2_ok=m2[m2[,2] %in% tmp1[,1],]

















