#' @title Violinplot_Shine
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param ppt_width default is 9
#' @param levels order of group
#' @param font_size font_size
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
Violinplot_Shine <- function(ppt_width = 9,
                          font_size = 20,
                          levels = c("Normal","CHD")

){
  pacman::p_load(ggplot2,ggsci)
  data <- read.csv("data for boxplot.csv")
  data$group <- factor(data$group, levels = levels)
  data$metabolites <- fct_inorder(data$metabolites)
  s <- ggplot(data,aes(x=metabolites,y = abundance))+
    geom_violin(aes(fill = group))+
    theme_bw()+
    scale_fill_npg()+
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1,size=font_size),
          legend.title=element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size = font_size),
          axis.text.y = element_text(size = font_size),
          axis.title.x = element_text(size = 0),#the font size of axis title
          axis.title.y = element_text(size = font_size)
    )
  save(s,file = 'violinplot.Rda')
  ggsave('boxplot.pdf',height=7,width=9)
  #export::graph2ppt(x=s,file='violinplot.pptx',height=12,width=ppt_width)
}




