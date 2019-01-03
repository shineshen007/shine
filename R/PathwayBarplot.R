#' pathway enrichment barplot
#' @export
#' @param row the rows be selected to plot,default is 30.
PathwayBarplot <- function(row = 30#draw the first 30 pathways
                           ){
  cat("Import data...\n")
  data <- read.csv("data.csv")
  data <- data[1:row,]
  group <- ifelse(data$p < 0.05,"sig", "not sig")
  pb<- ggplot2::ggplot(data,ggplot2::aes(reorder(pathway,-p),-log10(p)))+##-p control the order
    ggplot2::geom_bar(aes(fill=group),stat = "identity",position="dodge",width=0.8)+
    ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
          panel.background = element_blank(),axis.line = element_line(colour = "black"),
          legend.position = "none",axis.text.y = element_text(size = 14),#the font size of axis
          axis.title.x = element_text(size = 14),#the font size of axis title
          axis.title.y = element_text(size = 14))+
    ggplot2::geom_hline(yintercept= 1.30103,#draw a line at p=0.05
               lty=4,col="grey21",lwd=1)+
    ggplot2::labs(title="Pathway Enrichment Analysis")+
    ggplot2::coord_flip()+
    ggplot2::xlab('Pathway')+
    ggplot2::ggsave("PathwayBarplot.tiff", width = 12, height = 8)

}





