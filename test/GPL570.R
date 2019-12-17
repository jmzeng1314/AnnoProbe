library(AnnoProbe)
ids=idmap('GPL570',type = 'soft')
head(ids)
ids=ids[nchar(ids[,2])>1,]
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
sort(table(ids$biotypes))
