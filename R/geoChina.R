##' Download expression dataset by GSE id
##'
##' \code{geoChina} will download the expression matrix and phenotype data as ExpressionSet format
##'      from cloud in mainland China,
##'      it's a alternative method for getGEO function from GEOquery package.
##'      geoChina('gse1009') is the same as eSet=getGEO('gse1009', getGPL = F)
##'
##' @param gse input GSE id, such as GSE1009, GSE2546, gse1009.
##' @param mirror "tencent" only for now.
##' @param destdir The destination directory for data downloads.
##' @return a list of ExpressionSet, which contains the  expression matrix and phenotype data
##' @importFrom utils download.file
##' @importClassesFrom Biobase ExpressionSet
##' @examples
##' \dontrun{
##' geoChina('GSE1009',destdir=tempdir())
##' }
##' @export geoChina

geoChina <- function(gse='GSE2546',mirror='tencent',destdir=getwd()){
  # eSet=getGEO('GSE2546', destdir=".", AnnotGPL = F, getGPL = F)
  # http://49.235.27.111/GEOmirror/GSE2nnn/GSE2546_eSet.Rdata
  # gse='GSE2546';mirror='tencent'
  gse=toupper(gse)
  if(!gse %in% series.accession){
    stop('Your GSE may not be expression by array, or even not a GSE')
  }
  down=ifelse(as.numeric(gsub('GSE','',gse))<1000,
              paste0('/GEOmirror/GSEnnn/',gse,
                     '_eSet.Rdata'),
              paste0('/GEOmirror/',
                     gsub('[0-9][0-9][0-9]$','nnn',gse),'/',gse,
                     '_eSet.Rdata'))

  if(mirror=='tencent'){
    up='http://49.235.27.111'
  }
  OS <- .Platform$OS.type

  if (OS == "unix"){
    tpf=paste0(destdir,"/", gse, '_eSet.Rdata') # MAC file path
  } else if (OS == "windows"){
    tpf=paste0(destdir,"\\", gse, '_eSet.Rdata') # windows file path
  } else {
    stop("ERROR: OS could not be identified")
  }

  download.file(paste0(up,down),tpf,mode = "wb")
  suppressWarnings(load(tpf))
  # getGEO('GSE2546', destdir=".", AnnotGPL = F, getGPL = F)
  message(paste0("file downloaded in ",destdir,'\nyou can also use getGEO from GEOquery, by \ngetGEO(',
               shQuote(gse),
               ', destdir=".", AnnotGPL = F, getGPL = F)'
               ))
  return(gset)
}
