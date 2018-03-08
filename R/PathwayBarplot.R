###pathway enrichment barplot
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




