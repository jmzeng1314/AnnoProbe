##' Annotate gene IDs according to GTF files in gencode
##'
##' annoGene will return a data.frame of gene information or write them to a file (csv or html format).
##' The user should set a list of genes to be annotated, with  "ENSEMBL" or "SYMBOL" style.
##'
##' @param IDs a list of genes
##' @param ID_type the type of input IDs, should be "ENSEMBL" or "SYMBOL"
##' @param species choose human or mouse, or rat, default: human
##' @param out_file the filename, should be ".csv" or ".html".
##' @importFrom DT datatable saveWidget
##' @importFrom methods hasArg
##' @importFrom utils write.csv
##' @return a dataframe which columns contain genesymbol, biotypes, ensembl ids and the positions of genes
##' @examples
##' IDs <- c("DDX11L1", "MIR6859-1", "OR4G4P", "OR4F5")
##' ID_type = "SYMBOL"
##' annoGene(IDs, ID_type)
##' \donttest{
##' annoGene(IDs, ID_type,out_file = tempfile(fileext = ".html"))
##' annoGene(IDs, ID_type,out_file = tempfile(fileext = ".csv"))
##' }

##' @export
annoGene <- function(IDs,ID_type,species='human',out_file){
  if(length(unique(IDs))<1){
    stop("You should give me some genes to be annotated!!!")
  }
  if(!ID_type %in% c("ENSEMBL" ,"SYMBOL")){
    stop("We only accept  ENSEMBL or SYMBOL !!!")
  }
  if(species=='human'){
    GTF <- humanGTF
  }else if(species=='mouse'){
    GTF <- mouseGTF
  }else if(species=='rat'){
    GTF <- ratGTF
  }else{
    stop("We only accept human or mouse, or rat, ")
  }
  res <- GTF[eval(parse(text=paste0("GTF$",ID_type))) %in% IDs, ]

  missIds <- IDs[!(IDs %in% eval(parse(text=paste0("res$",ID_type))))]
  missIdsPercentage = round((length(missIds)/length(IDs))*100,2)
  if(length(missIds)!=0){
    warning(
      paste0(missIdsPercentage ,"% of input IDs are fail to annotate... ")
      # example: 5.29% of input gene IDs are fail to map...
    )
  }
  if (hasArg(out_file)) {
    results=res
    if(grepl('.html$',out_file)){
      Ensembl_prefix <- "https://asia.ensembl.org/Homo_sapiens/Gene/Summary?g="
      href = paste0(Ensembl_prefix, results$ENSEMBL)
      results$ENSEMBL = paste0("<b><a target=\"_black\" href=", shQuote(href), ">", results$ENSEMBL, "</a></b>")

      symbol_prefix <- "http://www.ncbi.nlm.nih.gov/gene?term="
      href = paste0(symbol_prefix, results$SYMBOL)
      results$SYMBOL = paste0("<b><a target=\"_black\" href=", shQuote(href), ">", results$SYMBOL, "</a></b>")

      y <- DT::datatable(results, escape = F, rownames = F)
      DT::saveWidget(y,file = out_file)
    }else if(grepl('.csv$',out_file)){
      write.csv(results,file =out_file )
    }else{
      stop("We only accept  csv or html format !!!")
    }

  }
  return(res)
}

