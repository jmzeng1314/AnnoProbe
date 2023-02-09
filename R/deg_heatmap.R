##' draw a heatmap for DEG result
##'
##' \code{deg_heatmap} will draw a heatmap for you.
##'
##' @param deg the result from limma.
##' @param genes_expr the expression matrix
##' @param  group_list, a vector
##' @param topn the number of genes in heatmap, default:20
##' @import ggplot2
##' @importFrom pheatmap pheatmap
##' @importFrom utils head tail
##' @return a ggplot2 style figure.
##' @examples
##' attach(GSE27533)
##' deg_heatmap(DEG,genes_expr,group_list)
##' @export
deg_heatmap <- function(deg,genes_expr,group_list,topn=20){
  x=deg[,1]
  names(x)=rownames(deg)
  cg=c(names(head(sort(x),topn)),
       names(tail(sort(x),topn)))
  n=t(scale(t(genes_expr[cg,])))
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  ac=data.frame(group_list=group_list)
  rownames(ac)=colnames(n)
  pheatmap(n,show_colnames =F,show_rownames = T,
           annotation_col=ac)
}
