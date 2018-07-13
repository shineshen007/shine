#' @title StaAnalysis
#' @description a function can generate PCA,PLSDA,heatmap,
#' s plot of foldchange and volcano plot, also can calculate vip value.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns.
#' @param sample.info a dataframe include sample.name,injection.order,
#' class,batch and group columns.
#' @param group group set.
#' @param p.cutoff default is 0.05.
#' @param splot default is FALSE.
#' @param pcorrect default is TRUE.
#' @param unitest t.test or wilcox.test.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' ##create a folder for Shine demo
#' dir.create("Demo for Shine")
#' setwd("Demo for Shine")
#' ##load the demo data
#' data(data, package = "Shine")
#' data(sample.info, package = "Shine")
#' ##export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.info, "sample.info.csv", row.names = FALSE)
#' ##run StaAnalysis
#' StaAnalysis(data = data,sample.info = sample.info,group = c("G","M"),
#' pcorrect = F)
#' }
StaAnalysis <- function(data = NULL,sample.info = NULL,p.cutoff = 0,
                       group = c("case","control"),
                       splot = FALSE,unitest =c("t.test","wilcox.test"),pcorrect = TRUE){
  cat("Analyzing data...\n")
  require(mixOmics);require(data.table)
  require(ggrepel);  require(pheatmap)

  cat("Import data...\n")
  data <- fread("data for sta.csv")
  data <- setDF(data)
  sample.info <- read.csv("sample.info.csv")


  ##create a folder for analysis
  path <-getwd()
  dir.create("StaAnalysis")
  setwd("StaAnalysis")
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  qc.name<-sample.info$sample.name[sample.info$class=="QC"]

  sample<-data[,match(sample.name,colnames(data))]
  qc<-data[,match(qc.name,colnames(data))]

  data_pfc<- as.matrix(cbind(qc,sample))

  name <- as.character(data[,"name"])

  class<- sample.info[,"group"]

  group1.index <- which(class == group[1])
  group2.index <- which(class == group[2])
  sample.index <- which(sample.info$class=="Subject")

  cat("Calculate Foldchange...\n")
  #must use the data_pfc,because the index include the qc in sampl.info
  fc <- apply(data_pfc,1,function(x) {
    median(x[group2.index]+0.1)/ median(x[group1.index]+0.1)
  })

  cat("Calculate P value...\n")

    test <- apply(data_pfc, 1, function(x) {
      unitest(x[group1.index], x[group2.index])
    })
    p <- unlist(lapply(test, function(x)
      x$p.value))

  if(pcorrect){
  p <- p.adjust(p = p, method = "fdr",n=length(p))
}
  cat("Draw PCA plot...\n")
  ###PCA
  png(file="PCA.png", width = 1200, height = 1000,res = 56*2)
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

  cat("Draw PLSDA plot...\n")
  ###PLS-DA
  png(file="PLSDA.png", width = 1200, height = 1000,res = 56*2)
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
            style="ggplot2",
            title = 'PLS-DA')
  dev.off()

  cat("Calculate VIP...\n")
  ###VIP
  vip<-vip(plsda.datatm)
  write.csv(vip,"VIP.csv",row.names = F)

  cat("Draw Volcano plot...\n")
  #### volcano plot(double line)###data only contain 3 columns:name,p,fc
  f<-as.data.frame(fc)
  pvalue<-as.data.frame(p)
  data_vol<-cbind(name,f,pvalue)
  write.csv(data_vol,"vol.csv",row.names = F)
  vol<-read.csv("vol.csv")
  fc<- vol$fc
  p<- vol$p
  group1<-group[1]
  group2<-group[2]
  Significant<- as.factor(ifelse(p < 0.05 & abs(log2(fc)) > 1,
                                 ifelse(log2(fc) < -1,
                                        "Down","Up"),"Not Sig"))
  png(file="volcano plot.png", width = 1200, height = 1000,res = 56*2)
  volc <- ggplot(vol, aes(x = log2(fc), y = -log10(p)))+
    geom_point(aes(color = Significant)) +
    scale_color_manual(values = c("green", "grey","red")) +
    annotate("text",x=xlim[2]-1,y=quantile(-log10(p)),label=group2)+
    annotate("text",x=xlim[2]-1.5,y=quantile(-log10(p)),label=group1)+
    theme_bw(base_size = 16) +
    geom_vline(xintercept=c(-1,1),
               lty=4,col="orange",lwd=1)+ # 在x轴-1.5与1.5的位置画两根竖线
    geom_hline(yintercept = -log10(0.05),
               lty=4,col="orange",lwd=1)+ #在p value 0.05的位置画一根横线
    labs(x="log2 (Fold change)",
         y="-log10 (p-value)",
         title="Volcano plot")+
    xlim(c(-5, 5))+
    geom_text_repel(
      data = subset(vol, p < p.cutoff&abs(log2(fc))>1),###fc的绝对值大于1
      max.iter = 100000,
      aes(label = name),
      size = 4,
      box.padding = 0.25,
      point.padding = 0.3
    )
  plot(volc)
  dev.off()

  if(splot){
  cat("Draw S plot of foldchange...\n")
  png(file="S plot of foldchange.png", width = 900, height = 800,res = 56*2)
  ##S plot of foldchange
  index<-c(1:nrow(f))
  datas <- cbind(index,f)
  Significant_s<- as.factor(ifelse(abs(log2(datas$fc)) > 0.26,
                                 ifelse(log2(datas$fc) < -0.27,
                                        "Down","Up"),"Not Sig"))
  splot <- ggplot(datas, aes(x = reorder(index,fc), y = log2(fc)))+
    geom_point(aes(color = Significant_s)) +
    scale_color_manual(values = c("green", "grey","red"))+
    labs(title="S plot of foldchange")+
    xlab('Index')
  plot(splot)

  dev.off()
}
  ##back origin work directory
  setwd(path)

  cat("StaAnalysis is done\n")

}
