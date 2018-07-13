#' @title volcano
#' @description a function can generate volcano plot.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns.
#' @param sample.info a dataframe include sample.name,injection.order,
#' class,batch and group columns.
#' @param group group set.
#' @param p.cutoff default is 0.05.
#' @param pcorrect default is TRUE.
#' @param doubleline default is TRUE.
#' @param singleline default is FALSE.
#' @param unitest t.test or wilcox.test.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
volcano <- function(data = NULL,sample.info = NULL,p.cutoff = 0,
                       group = c("case","control"),pcorrect = TRUE,
                    singleline = TRUE,xlim=c(-5,5),
                   doubleline = FALSE,unitest =c("t.test","wilcox.test")){
  require(data.table)
  cat("Import data...\n")
  data <- fread("data.csv")
  data <- setDF(data)
  sample.info <- read.csv("sample.info.csv")
  cat("Analyzing data...\n")
  require(mixOmics)
  require(ggrepel);  require(gplots)
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
  if(doubleline){
  Significant<- as.factor(ifelse(p < 0.05 & abs(log2(fc)) > 1,
                                 ifelse(log2(fc) < -1,
                                        "Down","Up"),"Not Sig"))
  png(file="volcano plot doubleline.png", width = 1200, height = 1000,res = 56*2)
  volc1 <- ggplot(vol, aes(x = log2(fc), y = -log10(p)))+
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
    xlim(xlim)+
    geom_text_repel(
      data = subset(vol, p < p.cutoff&abs(log2(fc))>1),###fc的绝对值大于1
      max.iter = 100000,
      aes(label = name),
      size = 4,
      box.padding = 0.25,
      point.padding = 0.3
    )
  plot(volc1)
  dev.off()
  }
  if(singleline){
    Significant<- as.factor(ifelse(p < 0.05 & abs(log2(fc)) > 0,
                                   ifelse(log2(fc) < 0,
                                          "Down","Up"),"Not Sig"))
    png(file="volcano plot sinleline.png", width = 1200, height = 1000,res = 56*2)
    volc2 <- ggplot(vol, aes(x = log2(fc), y = -log10(p)))+
      geom_point(aes(color = Significant)) +
      scale_color_manual(values = c("green", "grey","red")) +
      annotate("text",x=xlim[2]-1,y=quantile(-log10(p)),label=group2)+
      annotate("text",x=xlim[2]-1.5,y=quantile(-log10(p)),label=group1)+
      theme_bw(base_size = 16) +
      geom_vline(xintercept= 0,
                 lty=4,col="orange",lwd=1)+ # 在x轴-1.5与1.5的位置画两根竖线
      geom_hline(yintercept = -log10(0.05),
                 lty=4,col="orange",lwd=1)+ #在p value 0.05的位置画一根横线
      labs(x="log2 (Fold change)",
           y="-log10 (p-value)",
           title="Volcano plot")+
      xlim(xlim)+
      geom_text_repel(
        data = subset(vol, p < p.cutoff),###fc的绝对值大于1
        max.iter = 100000,
        aes(label = name),
        size = 4,
        box.padding = 0.25,
        point.padding = 0.3
      )
    plot(volc2)
    dev.off()
  }

}
