rm(list=ls())
options(stringsAsFactors = F)
load('../idmap1/data/p2s_df.rda')
library(usethis)
head(p2s_df)
ns=lapply(split(p2s_df,p2s_df$gpl), function(x){
  # x=tmp[[1]]
  ids=x[,-3]
  p=paste0(x[1,3],'_bioc')
  assign(p,ids, envir = .GlobalEnv)
  return(as.name(p))
})
ns
library(usethis)
do.call("use_data",ns)


rm(list=ls())
options(stringsAsFactors = F)
load('~/Documents/GPL/p2s_list_from_soft.Rdata')
length(p2s_list)
names(p2s_list)
head(p2s_list[[1]])
ns=lapply(p2s_list, function(x){
  # x=tmp[[1]]
  ids=x[,-3]
  assign(x[1,3],ids, envir = .GlobalEnv)
  return(as.name(x[1,3]))
})
ns
library(usethis)
do.call("use_data",ns)
# ls GPL*|while read id ;do ( cp  $id ${id%%.*}_soft.rda);done
# scp  *_soft.rda  jmzeng@49.235.27.111:/project/GEOmirror/GPL

## 把hisat2比对和bedtools依据gtf注释后的每个gpl的探针对应基因load进去
## 读取文件夹下面，全部的.probe2gene结尾的文件。
options(stringsAsFactors = F)
fs=list.files(pattern = 'probe2gene')
fs
p2s_list_pipe <- lapply(fs, function(x){
  a=read.table(file.path('./',x))
  colnames(a)=c('probe_id','symbol')
  a$gpl=paste0(strsplit(x,'_')[[1]][1],'_pipe')
  return(a)
})
lapply(p2s_list_pipe, head)

if(F){
  library(RSQLite)
  sqlite    <- dbDriver("SQLite")
  con <- dbConnect(sqlite,"probes_pipeline.sqlite") # makes a new file
  lapply(p2s_list_pipe, function(x){
    ids=x[,-3]
    print(x[1,3])
    dbWriteTable(con,x[1,3],ids,row.name=F,overwrite=T)
  })
  # it's too big, more than 200 Mb, so quit.
}

rm(list=ls())
options(stringsAsFactors = F)
load('p2s_list_pipe_lncRNA.Rdata')
ns=lapply(p2s_list_pipe, function(x){
  # x=tmp[[1]]
  ids=x[,-3]
  assign(x[1,3],ids, envir = .GlobalEnv)
  return(as.name(x[1,3]))
})
ns
library(usethis)
do.call("use_data",ns )


rm(list=ls())
options(stringsAsFactors = F)
gpl_list=read.csv('gpl_list.csv')
usethis::use_data(gpl_list)
# scp  *_soft.rda  jmzeng@49.235.27.111:/project/GEOmirror/GPL

rm(list=ls())
options(stringsAsFactors = F)
gpl_list=read.table('test/exists_anno.txt',header = F,sep='\t')
usethis::use_data(gpl_list,overwrite = T)
exists_anno_list=gpl_list[,1]
exists_anno_list=unique(exists_anno_list)
usethis::use_data(exists_anno_list,overwrite = T)




