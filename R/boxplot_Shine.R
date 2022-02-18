#' @title Boxplot_Shine
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param levels order of group
#' @param font_size font_size
#' @param height 7
#' @param width 9
#' @return  All the results can be got form other functions and instruction.
#' @export
Boxplot_Shine <- function(font_size = 20,
                          levels =c ("NU","HU","Gout"),
                          height=7,width=9

){
  pacman::p_load(ggpubr,ggsci)
  data <- read.csv("data for boxplot.csv")
  data$group <- factor(data$group, levels = levels)
  data$metabolites <- fct_inorder(data$metabolites)
  s <- ggpubr::ggboxplot(data = data,x='metabolites',y = 'abundance',color = "group",
                         add = "jitter"
                         #,ylab = ylab
  )+scale_color_nejm()+
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1,size=font_size),
          legend.title=element_blank(),
          legend.text = element_text(size = font_size),
          axis.text.y = element_text(size = font_size),
          axis.title.x = element_text(size = 0),#the font size of axis title
          axis.title.y = element_text(size = font_size)
    )+
    stat_compare_means(aes(group = group),label = "p.signif", label.x = 1.5)+
    ggsave('boxplot.pdf',height=height,width=height)
  save(s,file = 'boxplot.Rda')
  #export::graph2ppt(x=s,file='boxplot.pptx',height=7,width=ppt_width)
}




