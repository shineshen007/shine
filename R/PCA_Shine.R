#' @title PCA_Shine
#' @description a function can generate PCA.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param ind default is FALSE.
#' @param ellipse default is FALSE.
#' @param both default is FALSE.
#' @param neither default is TRUE.
#' @param QC draw the pca plot of QC.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
PCA_Shine <- function(ind = FALSE,ellipse = FALSE,
               both = FALSE,neither = TRUE,QC = FALSE){

  cat("Import data...\n")
  data <- data.table::fread("data.csv")
  data <- data.table::setDF(data)
  sample.info <- read.csv("sample.info.csv")
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  qc.name<-sample.info$sample.name[sample.info$class=="QC"]

  sample<-data[,match(sample.name,colnames(data))]
  qc<-data[,match(qc.name,colnames(data))]

  data_pfc<- as.matrix(cbind(qc,sample))

  name <- as.character(data[,"name"])

  class<- as.factor(sample.info[,"group"])
  batch<- as.factor(sample.info$batch)
  colour<-c("Turquoise3","Gold1","Firebrick1","Purple4","NavyBlue","LightSlateBlue",
            "Magenta","HotPink2","DeepSkyBlue2","Green1")#10 color
  lev<-length(levels(class))
  levb<-length(levels(batch))
  cl<-colour[1:lev]
  clo<-colour[1:levb]
  shape <-c(15:18,0,1,2,5,6)
  pch<-shape[1:lev]
  pchb<-shape[1:levb]

  if(QC){
    cat("Draw PCA plot of QC...\n")
    ###PCA"
    n <- which(sample.info$group=="QC")
    #png(file="PCA qc both.png", width = 1200, height = 1000,res = 56*2)
    temp<-qc
    pca<- mixOmics::pca(t(temp), ncomp=2, scale=T)
    pcap<- mixOmics::plotIndiv(pca,
                    group = as.factor(sample.info$batch[n]),
                    ind.names = T,###label
                    ellipse = T,###confidence interval
                    legend =TRUE,
                    cex=1.6,
                    style="graphics",
                    abline = T,
                    title = 'PCA')
    export::graph2ppt(file='PCA qc.png',height=7,width=9)
    #dev.off()

    #png(file="PCA qc neither.png", width = 1200, height = 1000,res = 56*2)
    temp<-qc
    pca<- mixOmics::pca(t(temp), ncomp=2, scale=T)
    pcap<-mixOmics::plotIndiv(pca,
                    group = as.factor(sample.info$batch[n]),
                    ind.names = F,###label
                    ellipse = F,###confidence interval
                    legend =TRUE,
                    pch = pchb,
                    cex=1.6,
                    style="graphics",
                    abline = T,
                    title = 'PCA')
    export::graph2ppt(file='PCA qc.png',height=7,width=9,append=TRUE)

    #dev.off()
  }

  cat("Draw PCA plot...\n")
  if(ind){
  ###PCA"
  #png(file="PCA ind.png", width = 1200, height = 1000,res = 56*2)
  temp<-data_pfc
  pca<-mixOmics::pca(t(temp), ncomp=2, scale=T)
  pcap<-mixOmics::plotIndiv(pca,
            group = sample.info$group,
            ind.names = T,###label
            ellipse = F,###confidence interval
            legend =TRUE,
            cex=1.6,
            style="graphics",
            abline = T,
            title = 'PCA')
  export::graph2ppt(file='PCA ind.png',height=7,width=9)
  #dev.off()
  }
  if(ellipse){
    ###PCA
    #png(file="PCA ellipse.png", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<-mixOmics::pca(t(temp), ncomp=2, scale=T)
    pcap<-mixOmics::plotIndiv(pca,
                    group = sample.info$group,
                    ind.names = F,###label
                    ellipse = T,###confidence interval
                    legend =TRUE,
                    col=cl,
                    pch = pch,
                    point.lwd=3,
                    cex=1.6,
                    style="graphics",
                    abline = T,
                    title = 'PCA')
    export::graph2ppt(file='PCA ellipse.png',height=7,width=9)
    #dev.off()
  }
  if(both){
    ###PCA
    #png(file="PCA both.png", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<-mixOmics::pca(t(temp), ncomp=2, scale=T)
    pcap<-mixOmics::plotIndiv(pca,
                    group = sample.info$group,
                    ind.names = T,###label
                    ellipse = T,###confidence interval
                    legend =TRUE,
                    style="graphics",
                    cex=1.6,
                    abline = T,
                    title = 'PCA')
    export::graph2ppt(file='PCA both.png',height=7,width=9)
    #dev.off()
  }
  if(neither){
    ###PCA
    #png(file="PCA neither.png", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<-mixOmics::pca(t(temp), ncomp=2, scale=T)
    pcap<-mixOmics::plotIndiv(pca,
                    group = sample.info$group,
                    ind.names = F,###label
                    ellipse = F,###confidence interval
                    legend =TRUE,
                    col=cl,
                    point.lwd=3,
                    pch = pch,
                    cex=1.6,
                    style="graphics",
                    abline = T,
                    title = 'PCA')
    export::graph2ppt(file='PCA neither.png',height=7,width=9)

    #dev.off()

  }
  dev.off()
}
