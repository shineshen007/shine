#' @title PathwayAnalysis
#' @description a function to do PathwayAnalysis
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data data
#' @param nrow 5
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
PathwayAnalysis <- function(data='data pathway.csv',
                            nrow=5){
  rawdata <- read.csv(data)
  dir.create('pathway')
  setwd('pathway')
  #library(MetaboAnalystR)
  #library(ggplot2)
  mSet<-InitDataObjects("conc", "pathora", FALSE)
  cmpd.vec <- unique(unlist(strsplit(as.character(rawdata$compound.name),";")))

  mSet<-Setup.MapData(mSet, cmpd.vec);
  Sys.setlocale(category="LC_ALL", locale = "English_United States.1252")
  mSet<-CrossReferencing(mSet, "name");
  mSet<-CreateMappingResultTable(mSet)
  mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
  mSet<-SetMetabolomeFilter(mSet, F);
  mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
  mSet<-PlotPathSummary(mSet, "path_view_0_", "png", 72, width=NA)
  mSet<-SaveTransformedData(mSet)
  #
  data <- read.csv("pathway_results.csv")
  colnames(data)[5] <- 'p'
  row=sum(data$p<0.05)+nrow
  data <- data[1:row,]
  group <- ifelse(data$p < 0.05,"sig", "not sig")
  pb <- ggplot2::ggplot(data,ggplot2::aes(reorder(X,-p),-log10(p)))+##-p control the order
    ggplot2::geom_bar(aes(fill=group),stat = "identity",position="dodge",width=0.8)+
    scale_fill_manual(values = c('Turquoise3','Firebrick1'))+
    ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
                   panel.background = element_blank(),axis.line = element_line(colour = "black"),
                   legend.position = "none",axis.text.y = element_text(size = 14),
                   axis.text.x = element_text(size = 14),#the font size of axis
                   axis.title.x = element_text(size = 14),#the font size of axis title
                   axis.title.y = element_text(size = 14))+
    ggplot2::geom_hline(yintercept= 1.30103,#draw a line at p=0.05
                        lty=4,col="grey21",lwd=1)+
    ggplot2::labs(title="Pathway Enrichment Analysis")+
    ggplot2::coord_flip()+
    ggplot2::xlab('Pathway')
    save(pb,file = 'barplot.Rda')

    export::graph2ppt(x=pb,file='pathway.pptx',height=7,width=9)
    #####
    vol<-read.csv("pathway_results.csv") #%>%

    Impact<- vol$Impact
    p<- vol$Raw.p


    sqx <- sqrt(Impact);
    min.x<- min(sqx, na.rm = TRUE);
    max.x <- max(sqx, na.rm = TRUE);

    if(min.x == max.x){ # only 1 value
      max.x = 1.5*max.x;
      min.x = 0.5*min.x;
    }

    maxR <- (max.x - min.x)/40;
    minR <- (max.x - min.x)/160;
    radi.vec <- minR+(maxR-minR)*(sqx-min.x)/(max.x-min.x)
    # set background color according to y
    bg.vec <- heat.colors(length(p));
    op <- par(mar=c(6,5,2,3));
    plot(Impact, -log(p), type="n", axes=F, xlab="Pathway Impact", ylab="-log(p)");
    axis(1);
    axis(2);
    grid(col="blue");
    symbols(Impact, -log(p), add = TRUE, inches = F, circles = radi.vec, bg = bg.vec, xpd=T);
    export::graph2ppt(file='pathway.pptx',height=7,width=9,append=T)

    ##
    ph <- ggplot2::ggplot(vol, aes(x = Impact, y = -log(Raw.p)))+
      geom_point(size=radi.vec*300,color=bg.vec) +
      #scale_color_manual(values = colour)
      geom_vline(xintercept=0.75,
                 lty=4,col="black",lwd=0.5)+ # add vetical line
      geom_hline(yintercept = -log(0.05),
                 lty=4,col="black",lwd=0.5)+ #add hori line
      labs(x="Pathway Impact",
           y="-log10(p-value)")+
      ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
                     panel.background = element_blank(),axis.line = element_line(colour = "black"),
                     legend.position = "none",axis.text.y = element_text(size = 14),
                     axis.text.x = element_text(size = 14),#the font size of axis
                     axis.title.x = element_text(size = 14),#the font size of axis title
                     axis.title.y = element_text(size = 14))+
      ggrepel::geom_text_repel(
        data = subset(vol,Raw.p< 0.05),
        #max.iter = 10,
        aes(label = X),
        size = 4,
        box.padding = 0.25,
        point.padding = 0.3)
    export::graph2ppt(x=ph,file='pathway.pptx',height=7,width=9,append=T)
    setwd('..//')
}

