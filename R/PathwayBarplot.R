#' pathway enrichment barplot
#' @export
#' @param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns.
PathwayBarplot <- function(data){
  require(ggplot2)
  cat("Import data...\n")
  data <- read.csv("data.csv")
  data <- data[1:30,]
  group <- ifelse(data$p < 0.05,"sig", "not sig")
  pb<- ggplot(data,aes(reorder(pathway,-p),-log10(p)))+##-p control the order
    geom_bar(aes(fill=group),stat = "identity",position="dodge",width=0.8)+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
          panel.background = element_blank(),axis.line = element_line(colour = "black"),
          legend.position = "none",axis.text.y = element_text(size = 14),#the font size of axis
          axis.title.x = element_text(size = 14),#the font size of axis title
          axis.title.y = element_text(size = 14))+
    labs(title="Pathway Enrichment Analysis")+
    coord_flip()+
    xlab('Pathway')+
    ggsave("PathwayBarplot.tiff", width = 12, height = 8)

}





