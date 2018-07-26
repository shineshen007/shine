#' @title PCA_Shine
#' @description a function can generate PCA.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns.
#' @param sample.info a dataframe include sample.name,injection.order,
#' class,batch and group columns.
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
PCA_Shine <- function(data = NULL,sample.info = NULL,ind = FALSE,ellipse = FALSE,
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

  cat("Draw PCA plot...\n")
  if(ind){
  ###PCA
  png(file="PCA ind.png", width = 1200, height = 1000,res = 56*2)
  temp<-data_pfc
  pca<-pca(t(temp), ncomp=2, scale=T)
  pcap<-plotIndiv(pca,
            group = sample.info$group,
            ind.names = T,###label
            ellipse = F,###confidence interval
            legend =TRUE,
            style="graphics",
            abline = T,
            title = 'PCA')

  dev.off()
  }
  if(ellipse){
    ###PCA
    png(file="PCA ellipse.png", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<-pca(t(temp), ncomp=2, scale=T)
    pcap<-plotIndiv(pca,
                    group = sample.info$group,
                    ind.names = F,###label
                    ellipse = T,###confidence interval
                    legend =TRUE,
                    style="graphics",
                    abline = T,
                    title = 'PCA')

    dev.off()
  }
  if(both){
    ###PCA
    png(file="PCA both.png", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<-pca(t(temp), ncomp=2, scale=T)
    pcap<-plotIndiv(pca,
                    group = sample.info$group,
                    ind.names = T,###label
                    ellipse = T,###confidence interval
                    legend =TRUE,
                    style="graphics",
                    abline = T,
                    title = 'PCA')

    dev.off()
  }
  if(neither){
    ###PCA
    png(file="PCA neither.png", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<-pca(t(temp), ncomp=2, scale=T)
    pcap<-plotIndiv(pca,
                    group = sample.info$group,
                    ind.names = F,###label
                    ellipse = F,###confidence interval
                    legend =TRUE,
                    style="graphics",
                    abline = T,
                    title = 'PCA')

    dev.off()
  }
}
