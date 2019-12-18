rm(list = ls())
library(AnnoProbe)
ids=idmap('GPL570',type = 'soft')
head(ids)

table(nchar(ids[,2])>1)
# 虽然从GPL的soft文件里面可以看到全部的 54675 个探针，但其实是有 8894 无法注释到基因。
ids=ids[nchar(ids[,2])>1,]
# 注释到多个基因的只有 2796 个探针
ids1=ids[grepl('///',ids[,2]),]
ids2=ids[!grepl('///',ids[,2]),]

# 我觉得下面的函数写的很差，运行太慢
tmp = do.call(rbind,apply(ids1,1,function(x){
  x[1];x[2]
  data.frame(ID=x[1],symbol=strsplit(x[2],' /// ')[[1]])
})
)
ids=rbind(ids2,tmp)
anno=annoGene(ids$symbol,"SYMBOL")
ids=merge(ids,anno,by.x = 'symbol',by.y='SYMBOL',all.x = T)

tail(sort(table(ids$biotypes)))

ids=idmap('GPL570',type = 'bioc')
anno=annoGene(ids$symbol,"SYMBOL")
ids=merge(ids,anno,by.x = 'symbol',by.y='SYMBOL',all.x = T)
tail(sort(table(ids$biotypes)))

