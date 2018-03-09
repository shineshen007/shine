#' pathway enrichment barplot
#' @export
#' @param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns.
PathwayBarplot<- function(data = NULL){
  png(file="PathwayBarplot.png", width = 900, height = 800)
  pb<- ggplot(data,aes(reorder(pathway,-p),-log10(p)))+##-p control the order
    geom_bar(aes(colour=p<0.05),stat = "identity",position="dodge",width=0.8)+
    labs(title="pathway enrichment analysis")+
    coord_flip()+
    xlab('pathway')
  plot(pb)
  dev.off()
}




