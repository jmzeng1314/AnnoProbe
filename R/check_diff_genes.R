#' @title Check a list of genes how they show difference.
#'
#' @description How does a gene or a list of genes show difference between two group.
#' The boxplot or heatmap will be drawed.
#' just a wrap function of ggpubr and pheatmap.
#' @param gene A vector contains all gene ids of interest. Gene ids should
#' be gene symbol.
#' @param genes_expr An expression matrix, the rownames should be  gene symbol.
#' @param group_list A vector contains the group information of each samples in  expression matrix
#' @export
#' @importFrom ggpubr ggboxplot
#' @importFrom pheatmap pheatmap
#' @return A figure : boxplot or heatmap
#' @examples
#' attach(GSE95166)
#' check_diff_genes('LRCH3',genes_expr,group_list )
#' \donttest{
#' x=DEG$logFC
#' names(x)=rownames(DEG)
#' cg=c(names(head(sort(x),100)),  names(tail(sort(x),100)))
#' check_diff_genes(cg,genes_expr,group_list )
#' }


check_diff_genes <- function(gene,genes_expr,group_list ){
  if(length(gene)==1){

    if(! gene %in% rownames(genes_expr)){
      stop(paste0(gene,' in not in your expression matrix'))
    }
    df=data.frame(value=as.numeric(genes_expr[gene,]),
                  group=group_list)
    ggpubr::ggboxplot(df, "group", "value",
              color = "group", palette =c("#00AFBB", "#E7B800"),
              add = "jitter", shape = "group")
  }else{
    cg=gene
    cg=cg[cg %in%  rownames(genes_expr) ]
    warning(paste0('Only ',length(cg),' in ',length(gene),' genes are in your expression matrix'))
    if(length(cg)<1){
      stop('None of the gene in your expression matrix')
    }
    n=t(scale(t(genes_expr[cg,])))
    n[n>2]=2
    n[n< -2]= -2
    n[1:4,1:4]
    ac=data.frame(group_list=group_list)
    rownames(ac)=colnames(n)
    pheatmap::pheatmap(n,show_colnames =F,show_rownames = F,
             annotation_col=ac)
  }

}
