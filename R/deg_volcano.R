##' draw a volcano for DEG result
##'
##' \code{deg_volcano} will draw a volcano for you.
##'
##' @param need_deg should be 3 columns : gene, logFC, p.value(or p.adjust
##' @param style you can try 1 or 2, default: 1
##' @param p_thred default:0.05
##' @param logFC_thred default:1
##' @importFrom ggplot2 ggplot aes geom_point theme_set theme_bw xlab ylab ggtitle theme element_text scale_colour_manual
##' @importFrom ggpubr ggscatter
##' @importFrom utils head
##' @export
##' @return a ggplot2 style figure.
##' @examples
##' deg=GSE27533$DEG
##' need_deg=data.frame(symbols=rownames(deg), logFC=deg$logFC, p=deg$P.Value)
##' deg_volcano(need_deg,2)
##' \donttest{
##' deg_volcano(need_deg,1)
##' }
deg_volcano <- function(need_deg,style=1,p_thred=0.05,logFC_thred=1){
  # need_deg should be 3 columns : gene, logFC, p.value(or p.adjust)
  colnames(need_deg)=c('gene','logFC','p')
  if(!(is.numeric(need_deg$logFC) & is.numeric(need_deg$p))){
    stop('we only need a data.frame which should be 3 columns : gene, logFC, p.value(or p.adjust)')
  }

  if(style==1){
    if(! logFC_thred){
      logFC_thred <- with(need_deg,mean(abs( logFC)) + 2*sd(abs( logFC)) )

    }
    # logFC_thred=1

    need_deg$change = as.factor(ifelse(need_deg$p < p_thred & abs(need_deg$logFC) > logFC_thred,
                                       ifelse(need_deg$logFC > logFC_thred ,'UP','DOWN'),'NOT')
    )
    this_tile <- paste0('Cutoff for logFC is ',round(logFC_thred,3),
                        '\nThe number of up gene is ',nrow(need_deg[need_deg$change =='UP',]) ,
                        '\nThe number of down gene is ',nrow(need_deg[need_deg$change =='DOWN',])
    )
    # message(this_tile)
    g = ggplot(data=need_deg,
               aes(x=logFC, y=-log10(p),
                   color=change)) +
      geom_point(alpha=0.4, size=1.75) +
      theme_set(theme_set(theme_bw(base_size=20)))+
      xlab("log2 fold change") + ylab("-log10 p-value") +
      ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
      scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
    return(g)
  }

  if(style==2){
      # p_thred=0.05;logFC_thred=1
      need_deg$g=ifelse(need_deg$p > p_thred,'stable',
                  ifelse( need_deg$logFC > logFC_thred,'up',
                          ifelse( need_deg$logFC < -logFC_thred,'down','stable') )
      )

      need_deg$p = -log10( need_deg$p)
      # message(table(need_deg$g))
      p=ggscatter(need_deg, x = "logFC", y = "p", color = "g",size = 0.5,
                label = "gene", repel = T,
                label.select =head(need_deg$gene),
                palette = c("#00AFBB", "#E7B800", "#FC4E07") )
      return(p)

    }

  ## TODO:
}

