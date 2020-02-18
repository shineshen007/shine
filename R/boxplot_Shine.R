#' @title Boxplot_Shine
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param ppt_width default is 9
#' @param ngroup number of grouo
#' @param palette colour selection
#' @param levels order of group
#' @param font_size font_size
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
Boxplot_Shine <- function(ppt_width = 9,
                          ngroup = 3,
                          font_size = 16,
                          levels =c ("Normal","Hyperuricemia","Gout"),
                          palette = c('BottleRocket1', 'BottleRocket2', 'Rushmore1', 'Royal1', 'Royal2',
                                      'Zissou1', 'Darjeeling1', 'Darjeeling2', 'Chevalier1' , 'FantasticFox1' ,
                                      'Moonrise1', 'Moonrise2', 'Moonrise3', 'Cavalcanti1', 'GrandBudapest1',
                                      'GrandBudapest2', 'IsleofDogs1', 'IsleofDogs2')
){
  data <- read.csv("data for boxplot.csv")
  data$group <- factor(data$group, levels = levels)
  s <- ggpubr::ggboxplot(data = data,x='metabolites',y = 'abundance',color = "group",
                         palette = wesanderson::wes_palette(n=ngroup, name=palette),
                         add = "jitter"
                         #,ylab = ylab
  )+
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          legend.title=element_blank(),
          legend.text = element_text(size = font_size),
          axis.text.y = element_text(size = font_size),
          axis.title.x = element_text(size = font_size),#the font size of axis title
          axis.title.y = element_text(size = font_size)
          )+
    stat_compare_means(aes(group = group),label = "p.signif", label.x = 1.5)
  save(s,file = 'boxplot.Rda')
  export::graph2ppt(x=s,file='boxplot.pptx',height=7,width=ppt_width)
}




