#' @title PLSDA_Shine
#' @description a function can generate PCA,PLSDA,heatmap,
#' s plot of foldchange and volcano plot, also can calculate vip value.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param ind default is FALSE.
#' @param ellipse default is FALSE.
#' @param both default is FALSE.
#' @param neither default is TRUE.
#' @param group group info.
#' @param multiclass multiclass plsda
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
PLSDA_Shine <- function(ind = FALSE,ellipse = FALSE,
                 both = FALSE,neither = TRUE,
                 group = c("case","control"),multiclass=FALSE){
  require(mixOmics)
  require(data.table)
  cat("Import data...\n")
  data <- fread("data.csv")
  data <- setDF(data)
  sample.info <- read.csv("sample.info.csv")
  class<- sample.info[,"group"]

  cat("Draw PLSDA plot...\n")
  if(multiclass){
    sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
    sample<-data[,match(sample.name,colnames(data))]
    sample.index <- which(sample.info$class=="Subject")
    ###PLS-DA
    png(file="PLSDA.png", width = 1200, height = 1000,res = 56*2)
    datat<-sample
    datatm<-as.matrix(datat)
    XXt<-t(datatm)
    group_pls<-as.data.frame(sample.info$group)
    YY<-group_pls[sample.index,]
    plsda.datatm <-plsda(XXt, YY, ncomp = 2)
    pls <- plotIndiv(plsda.datatm,
                     ind.names = F,
                     ellipse = T,
                     point.lwd=3,#point line size
                     legend =TRUE,
                     style="graphics",
                     title = 'PLS-DA')
    dev.off()
    cat("Calculate VIP...\n")
    ###VIP
    vip<-vip(plsda.datatm)
    write.csv(vip,"VIP.csv",row.names = F)
  }

  group1.index <- which(class == group[1])
  group2.index <- which(class == group[2])
  sample.info<-sample.info[c(group1.index,group2.index),]
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  sample<-data[,match(sample.name,colnames(data))]

  sample.index <- which(sample.info$class=="Subject")

  cat("Draw PLSDA plot...\n")
  if(ind){
  ###PLS-DA
  png(file="PLSDA ind.png", width = 1200, height = 1000,res = 56*2)
  datat<-sample
  datatm<-as.matrix(datat)
  XXt<-t(datatm)
  group_pls<-as.data.frame(sample.info$group)
  YY<-group_pls[sample.index,]
  plsda.datatm <-plsda(XXt, YY, ncomp = 2)
  pls <- plotIndiv(plsda.datatm,
            ind.names = T,
            ellipse = F,
            legend =TRUE,
            style="graphics",
            title = 'PLS-DA')
  dev.off()
  }
  if(ellipse){
    ###PLS-DA
    png(file="PLSDA ellipse.png", width = 1200, height = 1000,res = 56*2)
    datat<-sample
    datatm<-as.matrix(datat)
    XXt<-t(datatm)
    group_pls<-as.data.frame(sample.info$group)
    YY<-group_pls[sample.index,]
    plsda.datatm <-plsda(XXt, YY, ncomp = 2)
    pls <- plotIndiv(plsda.datatm,
                     ind.names = F,
                     ellipse = T,
                     point.lwd=3,#point line size
                     legend =TRUE,
                     style="graphics",
                     title = 'PLS-DA')
    dev.off()
  }
  if(both){
    ###PLS-DA
    png(file="PLSDA both.png", width = 1200, height = 1000,res = 56*2)
    datat<-sample
    datatm<-as.matrix(datat)
    XXt<-t(datatm)
    group_pls<-as.data.frame(sample.info$group)
    YY<-group_pls[sample.index,]
    plsda.datatm <-plsda(XXt, YY, ncomp = 2)
    pls <- plotIndiv(plsda.datatm,
                     ind.names = T,
                     ellipse = T,
                     legend =TRUE,
                     style="graphics",
                     title = 'PLS-DA')
    dev.off()
  }
  if(neither){
    ###PLS-DA
    png(file="PLSDA neither.png", width = 1200, height = 1000,res = 56*2)
    datat<-sample
    datatm<-as.matrix(datat)
    XXt<-t(datatm)
    group_pls<-as.data.frame(sample.info$group)
    YY<-group_pls[sample.index,]
    plsda.datatm <-plsda(XXt, YY, ncomp = 2)
    pls <- plotIndiv(plsda.datatm,
                     ind.names = F,
                     ellipse = F,
                     point.lwd=3,#point line size
                     legend =TRUE,
                     style="graphics",
                     title = 'PLS-DA')
    dev.off()
  }
  cat("Calculate VIP...\n")
  ###VIP
  vip<-vip(plsda.datatm)
  write.csv(vip,"VIP.csv",row.names = F)
}
