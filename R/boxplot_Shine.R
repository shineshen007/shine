#' @title Boxplot_Shine
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param ppt_width default is 9
#' @param ngroup number of grouo
#' @param palette colour selection
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
Boxplot_Shine <- function(ppt_width=9,
                          ngroup=3,
                          palette = c('BottleRocket1', 'BottleRocket2', 'Rushmore1', 'Royal1', 'Royal2',
                                      'Zissou1', 'Darjeeling1', 'Darjeeling2', 'Chevalier1' , 'FantasticFox1' ,
                                      'Moonrise1', 'Moonrise2', 'Moonrise3', 'Cavalcanti1', 'GrandBudapest1',
                                      'GrandBudapest2', 'IsleofDogs1', 'IsleofDogs2')
){
  data <- read.csv("data for boxplot.csv")
  s <- ggpubr::ggboxplot(data = data,x='metabolites',y = 'abundance',color = "group",
                         palette = wesanderson::wes_palette(n=ngroup, name=palette),
                         add = "jitter"
                         #,ylab = ylab
  )+
    theme(axis.text.x=element_text(family="myFont2",face="bold",angle=45, hjust=1, vjust=1),
          legend.title=element_blank(),
          legend.text = element_text(size = 16),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),#the font size of axis title
          axis.title.y = element_text(size = 18)
          )+
    stat_compare_means(aes(group = group),label = "p.signif", label.x = 1.5)

  plot(s)
  export::graph2ppt(file='boxplot.pptx',height=7,width=ppt_width)
}




