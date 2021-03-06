#' @title PLSDA_Shine
#' @description a function can generate PCA,PLSDA,heatmap,
#' s plot of foldchange and volcano plot, also can calculate vip value.
#' @author Xia Shen
#' \email{qq951633542@163.com}
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
#'  #the data must be named as 'data.csv' and 'sample.info.csv'
#' PLSDA()
#' }
PLSDA_Shine <- function(ind = FALSE,ellipse = FALSE,
                  both = FALSE,neither = TRUE,
                  group = c("case","control"),multiclass=FALSE){

  cat("Import data...\n")
  data <- data.table::fread("data.csv")
  data <- data.table::setDF(data)
  sample.info <- read.csv("sample.info.csv")
  class<- sample.info[,"group"]


  if(multiclass){
    cat("Draw PLSDA plot fot multiclass...\n")
    sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
    sample<-data[,match(sample.name,colnames(data))]
    sample.index <- which(sample.info$class=="Subject")
    ###PLS-DA
    png(filename="PLSDA.png", width = 1200, height = 1000,res = 56*2)
    datat<-sample
    datatm<-as.matrix(datat)
    XXt<-t(datatm)
    group_pls<-as.data.frame(sample.info$group)
    colour<-c("SeaGreen","goldenrod","Firebrick","DimGrey","NavyBlue","LightSlateBlue",
              "Magenta","HotPink2","DeepSkyBlue2","Green1")#10 color
    lev<-length(levels(class))
    cl<-colour[1:lev]
    shape <-c(15:18,7:14)
    pch<-shape[1:lev]
    YY<-group_pls[sample.index,]
    plsda.datatm <-mixOmics::plsda(XXt, YY, ncomp = 2)
    pls <- mixOmics::plotIndiv(plsda.datatm,
                               ind.names = F,
                               ellipse = T,
                               pch = pch,
                               cex=1.6,
                               point.lwd=3,#point line size
                               legend =TRUE,
                               style="graphics",
                               title = 'PLS-DA')
    dev.off()
    cat("Calculate VIP...\n")
    ###VIP
    vip<-mixOmics::vip(plsda.datatm)
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
    png(filename="PLSDA ind.png", width = 1200, height = 1000,res = 56*2)
    datat<-sample
    datatm<-as.matrix(datat)
    XXt<-t(datatm)
    group_pls<-as.data.frame(sample.info$group)
    YY<-group_pls[sample.index,]
    plsda.datatm <-mixOmics::plsda(XXt, YY, ncomp = 2)
    pls <- mixOmics::plotIndiv(plsda.datatm,
                               ind.names = T,
                               ellipse = F,
                               cex=1.6,
                               legend =TRUE,
                               style="graphics",
                               title = 'PLS-DA')
    dev.off()
  }
  if(ellipse){
    ###PLS-DA
    png(filename="PLSDA ellipse.png", width = 1200, height = 1000,res = 56*2)
    datat<-sample
    datatm<-as.matrix(datat)
    XXt<-t(datatm)
    group_pls<-as.data.frame(sample.info$group)
    YY<-group_pls[sample.index,]
    plsda.datatm <-mixOmics::plsda(XXt, YY, ncomp = 2)
    pls <- mixOmics::plotIndiv(plsda.datatm,
                               ind.names = F,
                               ellipse = T,
                               pch = 16,
                               cex=1.6,
                               point.lwd=3,#point line size
                               legend =TRUE,
                               style="graphics",
                               title = 'PLS-DA')
    dev.off()
  }
  if(both){
    ###PLS-DA
    png(filename="PLSDA both.png", width = 1200, height = 1000,res = 56*2)
    datat<-sample
    datatm<-as.matrix(datat)
    XXt<-t(datatm)
    group_pls<-as.data.frame(sample.info$group)
    YY<-group_pls[sample.index,]
    plsda.datatm <-mixOmics::plsda(XXt, YY, ncomp = 2)
    pls <- mixOmics::plotIndiv(plsda.datatm,
                               ind.names = T,
                               ellipse = T,
                               cex=1.6,
                               legend =TRUE,
                               style="graphics",
                               title = 'PLS-DA')
    dev.off()
  }
  if(neither){
    ###PLS-DA
    png(filename="PLSDA neither.png", width = 1200, height = 1000,res = 56*2)
    datat<-sample
    datatm<-as.matrix(datat)
    XXt<-t(datatm)
    group_pls<-as.data.frame(sample.info$group)
    YY<-group_pls[sample.index,]
    plsda.datatm <-mixOmics::plsda(XXt, YY, ncomp = 2)
    pls <- mixOmics::plotIndiv(plsda.datatm,
                               ind.names = F,
                               pch = 16,
                               cex=1.6,
                               ellipse = F,
                               point.lwd=3,#point line size
                               legend =TRUE,
                               style="graphics",
                               title = 'PLS-DA')
    export::graph2ppt(file='PLSDA neither.png',height=7,width=9)
    dev.off()

  }
  cat("Calculate VIP...\n")
  ###VIP
  vip<-mixOmics::vip(plsda.datatm)
  write.csv(vip,"VIP.csv",row.names = F)
}
