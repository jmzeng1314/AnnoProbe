##' Get Probe Annotation
##'
##' \code{idmap} returns probe annotations for input gpl
##' @param gpl GPL(GEO platform) number, eg: GPL570
##' @param type source of probe anntation stored, one of "pipe", "bioc", "soft", default:"pipe"
##' @param mirror "tencent" only for now
##' @param destdir The destination directory for data downloads.
##' @return probe annotaions
##' @importFrom utils download.file data
##' @examples
##' ids=idmap('GPL570',destdir=tempdir())
##' \donttest{
##' ids=idmap('GPL570',type='soft',destdir=tempdir())
##' ids=idmap('GPL18084',type='pipe',destdir=tempdir())
##' }
##' @export
##'
idmap <- function(gpl='GPL570',type='bioc',mirror='tencent',destdir=getwd()){
  gpl=toupper(gpl)
  gpl_anno=paste(gpl,c('bioc','soft','pipe'),sep='_')
  if(mirror=='tencent'){
    up='http://49.235.27.111'
  }
  if(!checkGPL(gpl)){
    stop("This platform is not in our list, please use our shinyAPP to custom annotate your probe sequences, or ask us to process and then update the R package!")
  }else{
    tryCatch("exists_anno_list")
    gpl_anno=gpl_anno [gpl_anno %in% exists_anno_list]

    # if( paste(gpl, type,sep='_')  %in% exists_anno_list){
     if(T){
      tpf=paste0( paste(gpl, type,sep='_'),'.rda')

      OS <- .Platform$OS.type

      if (OS == "unix"){
        dt = paste0( destdir,"/", tpf) # MAC file path
      } else if (OS == "windows"){
        dt = paste0( destdir,"\\", tpf) # windows file path
      } else {
        stop("ERROR: OS could not be identified")
      }


      down=paste0('/GEOmirror/GPL/',tpf)
      download.file(paste0(up,down),dt,mode = "wb")
      message(paste0("file downloaded in ",destdir))
      load(dt)
      return(get(paste(gpl, type,sep='_')))
    }else{
      stop('We have that platform, but just offer other type of annotaion.')
    }
  }

}


##' Check whether the input gpl in our platform list or not
##' @param GPL GPL(GEO platform) number, eg: GPL570
##' @return returns a boolean value
##' @examples
##' checkGPL('GPL570')
##' checkGPL('GPL15314')
##' checkGPL('GPL10558')
##' @export
checkGPL <- function(GPL=NULL){
  if(length(GPL)==0){
    stop("please input GPL number")
  }
  GPLList <- getGPLList()
  flag = (GPL %in% GPLList[,1])
  return(flag)
}

##' Print GPL information
##' @param GPL GPL(GEO platform) number, eg: GPL570
##' @return print detail information of the input GEO platform
##' @examples
##' printGPLInfo('GPL93')
##' @export
printGPLInfo <- function(GPL=NULL){
  if(length(GPL)!=0){
    flag=checkGPL(GPL)
    if(!flag){
      stop("This platform is not in our list, please use our shinyAPP to custom annotate your probe sequences, or ask us to process and then update the R package!")

    }
    tryCatch("gpl_list")
    gpl_list <- gpl_list[gpl_list[,1]==GPL,]
  }else{
    gpl_list <- t(getGPLList())
  }
  return(t(gpl_list))
}


##' Get all GPL list in our package
##' \code{getGPLList} returns all the GPL number checklist stored in package
##' @return a data.frame which contains the gpl and name of array.
##' @export
getGPLList <- function(){
  tryCatch("gpl_list")
  GPLList <- get("gpl_list")
  return(GPLList[,1:2])
}

