#' @title StaAnalysis
#' @description a function can generate PCA,PLSDA,heatmap,
#' s plot of foldchange and volcano plot, also can calculate vip value.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param group group set.
#' @param xlim xlim.
#' @param p.cutoff default is 0.05.
#' @param splot default is FALSE.
#' @param pcorrect default is TRUE.
#' @param heatmap heatmap
#' @param unitest t.test or wilcox.test.
#' @param paired paired test in t.test and wilcox.test,default is FALSE.
#' @param h the height of group index,default is 0.2
#' @param PCA draw PCA plot or not
#' @param colv the colour of volcano
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @import ggplot2
#' @examples
#' \donttest{
#' }
StaAnalysis <- function(p.cutoff = 0,
                        group = c("case","control"),
                        splot = FALSE,
                        colv = c("SpringGreen3","grey","Firebrick1"),
                        unitest =c("t.test","wilcox.test"),
                        pcorrect = TRUE,
                        xlim = c(-3,3),
                        paired = FALSE,
                        h=0.2,
                        PCA = FALSE,
                        heatmap = FALSE){
  cat("Analyzing data...\n")
  cat("Import data...\n")
  data <- data.table::fread("data for sta.csv")
  data <- data.table::setDF(data)
  sample.info <- read.csv("sample.info.csv")

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
    median(x[group1.index]+0.1)/ median(x[group2.index]+0.1)
  })

  cat("Calculate P value...\n")

  if(unitest == "t.test"){
    test <- apply(data_pfc, 1, function(x) {
      t.test(x[group1.index], x[group2.index],paired = paired)
    })
    p <- unlist(lapply(test, function(x)
      x$p.value))
  }
  if(unitest == "wilcox.test"){
    test <- apply(data_pfc, 1, function(x) {
      wilcox.test(x[group1.index], x[group2.index],paired = paired)
    })
    p <- unlist(lapply(test, function(x)
      x$p.value))
  }

  if(pcorrect){
    p <- p.adjust(p = p, method = "fdr",n=length(p))
  }
  ##create a folder for analysis
  FolderName <- paste("StaAnalysis",group[1],sep = " ")
  FolderName <- paste(FolderName,group[2],sep = "&")
  dir.create(FolderName)
  file.copy('sample.info.csv',FolderName)
  setwd(FolderName)

  #parameter decision
  unitest <- match.arg(unitest)
  group <- group
  correct <- as.logical(pcorrect)
  p.cutoff <- as.numeric(p.cutoff)
  paired <- as.logical(paired)

  ##save parameters
  StaAnalysis.parameters <- c(
    unitest,
    paste(group, collapse = ","),
    correct,
    paired,
    p.cutoff
  )
  StaAnalysis.parameters <- data.frame(c(
    "unitest",
    "group",
    "correct",
    "paired",
    "p.cutoff"
  ),
  StaAnalysis.parameters, stringsAsFactors = FALSE)

  StaAnalysis <- rbind(StaAnalysis.parameters,c("Version", "Shine"))
  colnames(StaAnalysis) <- c('parameter', 'value')
  write.csv(StaAnalysis,"StaAnalysis.parameters.csv",row.names = F)

  if(PCA){
    cat("Draw PCA plot...\n")
    #png(file="PCA.png", width = 1200, height = 1000,res = 56*2)
    temp<-data_pfc
    pca<- mixOmics::pca(t(temp), ncomp=2, scale=T)
    pcap<- mixOmics::plotIndiv(pca,
                               group = sample.info$group,
                               ind.names = F,###label
                               ellipse = F,###confidence interval
                               legend =TRUE,

                               point.lwd=3,

                               cex=1.6,
                               style="graphics",
                               abline = T,
                               title = 'PCA')
    export::graph2ppt(x=pcap,file='data.pptx',height=7,width=9)
    #dev.off()
  }


  cat("Draw PLSDA plot...\n")
  ###PLS-DA
  #png(file="PLSDA.png", width = 1200, height = 1000,res = 56*2)
  sample.info1<-sample.info[c(group1.index,group2.index),]
  ###data preparation
  sample.name1<-sample.info1$sample.name[sample.info1$class=="Subject"]
  sample1<-data[,match(sample.name1,colnames(data))]
  sample.index <- which(sample.info1$class=="Subject")
  #begin
  datat<-sample1
  datatm<-as.matrix(datat)
  XXt<-t(datatm)
  group_pls<-as.data.frame(sample.info1$group)
  YY<-group_pls[sample.index,]
  plsda.datatm <- mixOmics::plsda(XXt, YY, ncomp = 2)
  pls <- as.ggplot(function()mixOmics::plotIndiv(plsda.datatm,
                                                 ind.names = F,
                                                 ellipse = T,
                                                 pch = 16,#point shape
                                                 cex=1.6,#point size
                                                 point.lwd=3,# line size
                                                 legend =TRUE,
                                                 style="ggplot2",
                                                 title = 'PLS-DA'))
  save(pls,file = 'plsda_neither.Rda')
  export::graph2ppt(x=pls,file='data.pptx',height=7,width=9,append=TRUE)

  cat("Calculate VIP...\n")
  ###VIP
  vip<- mixOmics::vip(plsda.datatm)
  write.csv(vip,"VIP.csv",row.names = F)

  cat("Draw Volcano plot...\n")
  #### volcano plot(double line)###data only contain 3 columns:name,p,fc
  f<-as.data.frame(fc)
  pvalue<-as.data.frame(p)
  data_vol<-cbind(name,f,pvalue)
  write.csv(data_vol,"vol.csv",row.names = F)
  data_pfc_vip<-cbind(data_vol,vip,data)
  write.csv(data_pfc_vip,"data_pfc_vip.csv",row.names = F)

  #volcano plot
  DataFilter()
  vol<-read_csv("vol.csv")
  fc<- vol$fc
  p<- vol$p
  group1<-group[1]
  group2<-group[2]
  Significant<- as.factor(ifelse(p < 0.05 & abs(log2(fc)) > 0.41,
                                 ifelse(log2(fc) < -0.41,
                                        "Down","Up"),"Not Sig"))
  #png(file=, width = 1200, height = 1000,res = 56*2)
  volc <- ggplot2::ggplot(vol, ggplot2::aes(x = log2(fc), y = -log10(p)))+
    ggplot2::geom_point(aes(color = Significant),size=3) +
    ggplot2::scale_color_manual(values = colv) +
    ggplot2::annotate("text",x=xlim[2]-1,y=quantile(-log10(p),0.9999)-h,label=group2)+
    ggplot2::annotate("text",x=xlim[2]-1,y=quantile(-log10(p),0.9999),label=group1)+
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::geom_vline(xintercept=c(-0.41,0.41),
                        lty=4,col="orange",lwd=1)+ #
    ggplot2::geom_hline(yintercept = -log10(0.05),
                        lty=4,col="orange",lwd=1)+ #
    labs(x="log2 (Fold change)",
         y="-log10 (p-value)",
         title="Volcano plot")+
    xlim(xlim)+
    ggrepel::geom_text_repel(
      data = subset(vol, p < p.cutoff&abs(log2(fc))>0.41),###
      max.iter = 100000,
      aes(label = name),
      size = 4,
      box.padding = 0.25,
      point.padding = 0.3
    )
  save(volc,file = 'volcano.Rda')
  #ggplot2::ggsave("volcano plot.png", width = 12, height = 8)
  export::graph2ppt(x=volc,file='data.pptx',height=7,width=9,append = TRUE)

  if(heatmap){
    write.csv(sample.info,'sample.info.csv',row.names = F)
    Heatmap_Shine(group = group)
  }
  if(splot){
    cat("Draw S plot of foldchange...\n")
    png(file="S plot of foldchange.png", width = 900, height = 800,res = 56*2)
    ##S plot of foldchange
    index<-c(1:nrow(f))
    datas <- cbind(index,f)
    Significant_s<- as.factor(ifelse(abs(log2(datas$fc)) > 0.41,
                                     ifelse(log2(datas$fc) < -0.41,
                                            "Down","Up"),"Not Sig"))
    splot <- ggplot2::ggplot(datas, aes(x = reorder(index,fc), y = log2(fc)))+
      geom_point(aes(color = Significant_s)) +
      scale_color_manual(values = c("SpringGreen3", "grey","Firebrick1"))+
      labs(title="S plot of foldchange")+
      xlab('Index')
    plot(splot)

    dev.off()
  }
  ##back origin work directory
  setwd("..//")

  cat("StaAnalysis is done\n")

}
