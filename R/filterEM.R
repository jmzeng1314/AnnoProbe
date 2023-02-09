##' Filter expression matrix based on annotation
##'
##' \code{filterEM} will annotate the probes in expression matrix and remove the duplicated gene symbols.
##' because there will be many probes mapped to same genes, we will only keep the max value one.
##' @param probes_expr is an expression matrix which rownames are probes of probe2gene and each column is a sample
##' @param probe2gene  the first column is probes and the second column is corresponding gene symbols
##' @return a expression matrix which has been filtered duplicated gene symbols
##' @importFrom utils head
##' @importFrom stats na.omit median
##' @examples
##' attach(GSE95166)
##' # head(probes_expr)
##' # head(probe2gene)
##' genes_expr <- filterEM(probes_expr,probe2gene)
##' # head(genes_expr)
##' @export
filterEM <- function(probes_expr,probe2gene){
  colnames(probe2gene) <- c("probeid","symbol")
  probe2gene$probeid=as.character(probe2gene$probeid)
  probe2gene$symbol=trimws(probe2gene$symbol)
  # head(probe2gene)

  message(paste0('input expression matrix is ',nrow(probes_expr),' rows(genes or probes) and ',ncol(probes_expr),' columns(samples).\n'))
  message(paste0('input probe2gene is ',nrow(probe2gene),' rows(genes or probes)\n'))

  probe2gene=na.omit(probe2gene)
  # if one probe mapped to many genes, we will only keep one randomly.
  probe2gene=probe2gene[!duplicated(probe2gene$probeid),]
  # 这个地方是有问题的，随机挑选一个注释进行后续分析。
  probe2gene = probe2gene[probe2gene$probeid %in% rownames(probes_expr),]

  message(paste0('after remove NA or useless probes for probe2gene, ',nrow(probe2gene),' rows(genes or probes) left\n'))

  #probes_expr <- exprs(eSet);dim(probes_expr)
  probes_expr <- as.data.frame(probes_expr)
  message(paste0('There are ',
  sum(rownames(probes_expr) %in% probe2gene$probeid),
  ' of ',nrow(probes_expr),' probes can be annotated.\n'))

  probes_expr=probes_expr[as.character(probe2gene$probeid),]
  # probes_expr[1:4,1:4]

  probe2gene$median=apply(probes_expr,1,median)
  probe2gene=probe2gene[order(probe2gene$symbol,probe2gene$median,decreasing = T),]
  probe2gene=probe2gene[!duplicated(probe2gene$symbol),]


  genes_expr=probes_expr[as.character(probe2gene$probeid),]
  rownames(genes_expr)=probe2gene$symbol
  # genes_expr[1:4,1:4]
  message(paste0('output expression matrix is ',nrow(genes_expr),' rows(genes or probes) and ',ncol(genes_expr),' columns(samples).'))
  # probes_expr['AGAP6',]
  return(genes_expr)
}





