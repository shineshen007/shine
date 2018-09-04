#' @title PCA_Shine
#' @description a function can generate PCA.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param ind default is FALSE.
#' @param ellipse default is FALSE.
#' @param both default is FALSE.
#' @param neither default is TRUE.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
PCA_Shine <- function(ind = FALSE,ellipse = FALSE,
               both = FALSE,neither = TRUE){
  require(mixOmics)
  require(data.table)
  cat("Import data...\n")
  data <- fread("data.csv")
  data <- setDF(data)
  sample.info <- read.csv("sample.info.csv")
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  qc.name<-sample.info$sample.name[sample.info$class=="QC"]

  sample<-data[,match(sample.name,colnames(data))]
  qc<-data[,match(qc.name,colnames(data))]

  data_pfc<- as.matrix(cbind(qc,sample))

  name <- as.character(data[,"name"])

  class<- sample.info[,"group"]
  colour<-c("SeaGreen","goldenrod","Firebrick","DimGrey","NavyBlue","LightSlateBlue",
            "Magenta","HotPink2","DeepSkyBlue2","Green1")#10 color
  lev<-length(levels(class))
  cl<-colour[1:lev]
  shape <-c(15:18,7:14)
  pch<-shape[1:lev]

  cat("Draw PCA plot...\n")
  if(ind){
  ###PCA"
  tiff(file="PCA ind.tiff", width = 1200, height = 1000,res = 56*2)
  temp<-data_pfc
  pca<-pca(t(temp), ncomp=2, scale=T)
  pcap<-plotIndiv(pca,
            group = sample.info$group,
            ind.names = T,###label
            ellipse = F,###confidence interval
            legend =TRUE,
            pch = pch,
            cex=1.6,
            style="graphics",
            abline = T,
            title = 'PCA')

  dev.off()
  }
  if(ellipse){
    ###PCA
    tiff(file="PCA ellipse.tiff", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<-pca(t(temp), ncomp=2, scale=T)
    pcap<-plotIndiv(pca,
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

    dev.off()
  }
  if(both){
    ###PCA
    tiff(file="PCA both.tiff", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<-pca(t(temp), ncomp=2, scale=T)
    pcap<-plotIndiv(pca,
                    group = sample.info$group,
                    ind.names = T,###label
                    ellipse = T,###confidence interval
                    legend =TRUE,
                    style="graphics",
                    pch = pch,
                    cex=1.6,
                    abline = T,
                    title = 'PCA')

    dev.off()
  }
  if(neither){
    ###PCA
    tiff(file="PCA neither.tiff", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<-pca(t(temp), ncomp=2, scale=T)
    pcap<-plotIndiv(pca,
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

    dev.off()
  }
}
