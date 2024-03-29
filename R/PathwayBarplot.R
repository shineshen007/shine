#' @title PathwayBarplot
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param color color
#' @param font_size font_size
#' @return  All the results can be got form other functions and instruction.
#' @export
PathwayBarplot <- function(font_size=26,
                           color = c("#3399FF","#FF3333")
){
  cat("Import data...\n")
  data <- read.csv("pathway_results.csv")
  row=sum(data$q.value<0.05)+3
  data <- data[1:row,]
  group <- ifelse(data$q.value < 0.05,"sig", "not sig")
  pb<- ggplot2::ggplot(data,ggplot2::aes(reorder(X,-q.value),-log10(q.value)))+##-p control the order
    ggplot2::geom_bar(aes(fill=group),stat = "identity",position="dodge",width=0.8)+
    ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
                   panel.background = element_blank(),axis.line = element_line(colour = "black"),
                   title=element_text(size=font_size),
                   legend.position = "none",axis.text.y = element_text(size = font_size),#the font size of axis
                   axis.text.x = element_text(size = font_size),
                   axis.title.x = element_text(size = font_size),#the font size of axis title
                   axis.title.y = element_text(size = font_size))+
    ggplot2::geom_hline(yintercept= 1.30103,#draw a line at p=0.05
                        lty=4,col="grey21",lwd=1)+
    scale_fill_manual(values = color)+
    ggplot2::labs(title="Pathway Enrichment Analysis")+
    xlab(NULL)+
    ylab('-Log10(pvalue)')+
    ggplot2::coord_flip()#+
  #ggplot2::xlab('Pathway')#+
  ggplot2::ggsave("PathwayBarplot.pdf", width = 18, height = 20)
  save(pb,file = 'barplot.Rda')
}











