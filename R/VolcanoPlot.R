#' @title volcano
#' @description a function can generate volcano plot.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param group group set.
#' @param xlim xlim.
#' @param p.cutoff default is 0.05.
#' @param pcorrect default is TRUE.
#' @param doubleline default is TRUE.
#' @param singleline default is FALSE.
#' @param unitest t.test or wilcox.test.
#' @param fc.cutoff default is 1
#' @param colour the colour of points
#' @param paired paired test in t.test and wilcox.test,default is FALSE.
#' @param h the height of group index,default is 0.2
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
volcano <- function(p.cutoff = 0,
                    group = c("case","control"),
                    colour = c("SpringGreen3", "grey","Firebrick1"),
                    pcorrect = TRUE,
                    singleline = TRUE,
                    xlim = c(-3,3),
                    fc.cutoff = 0.41,
                    doubleline = FALSE,
                    unitest =c("t.test","wilcox.test"),
                    paired = FALSE,
                    h=0.2){

  cat("Import data...\n")
  data <- data.table::fread("data.csv")
  data <- data.table::setDF(data)
  sample.info <- read.csv("sample.info.csv")
  cat("Analyzing data...\n")

  #parameter decision
  unitest <- match.arg(unitest)
  group <- group
  correct <- as.logical(pcorrect)
  paired <- as.logical(paired)
  p.cutoff <- as.numeric(p.cutoff)
  fc.cutoff <- as.numeric(fc.cutoff)


  ##save parameters
  Volcano.parameters <- c(
    unitest,
    paste(group, collapse = ","),
    correct,
    paired,
    p.cutoff,
    fc.cutoff
  )
  Volcano.parameters <- data.frame(c(
    "unitest",
    "group",
    "correct",
    "paired",
    "p.cutoff",
    "fc.cutoff"
  ),
  Volcano.parameters, stringsAsFactors = FALSE)

  Volcano <- rbind(Volcano.parameters,c("Version", "0.0.984"))
  colnames(Volcano) <- c('parameter', 'value')
  write.csv(Volcano,"Volcano.parameters.csv",row.names = F)

  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  qc.name<-sample.info$sample.name[sample.info$class=="QC"]

  sample<-data[,match(sample.name,colnames(data))]
  qc<-data[,match(qc.name,colnames(data))]

  data_pfc<- as.matrix(cbind(qc,sample))

  name <- as.character(data[,"name"])

  class<- sample.info[,"group"]

  group1.index <- which(class == group[1])#case
  group2.index <- which(class == group[2])#control

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

  cat("Draw Volcano plot...\n")
  #### volcano plot(double line)###data only contain 3 columns:name,p,fc
  f<-as.data.frame(fc)
  pvalue<-as.data.frame(p)
  data_vol<-cbind(name,f,pvalue)
  write.csv(data_vol,"vol.csv",row.names = F)
  vol<-read.csv("vol.csv")
  fc<- vol$fc
  p<- vol$p
  group1<-group[1]#group index on the plot
  group2<-group[2]
  if(doubleline){
  Significant<- as.factor(ifelse(p < 0.05 & abs(log2(fc)) > 0.41,
                                 ifelse(log2(fc) < -0.41,
                                        "Down","Up"),"Not Sig"))
  png(file="volcano plot doubleline.png", width = 1200, height = 1000,res = 56*2)
  volc1 <- ggplot2::ggplot(vol, aes(x = log2(fc), y = -log10(p)))+
    geom_point(aes(color = Significant),size=3) +
    scale_color_manual(values = colour) +
    annotate("text",x=xlim[2]-1,y=quantile(-log10(p),0.9999)-h,label=group2)+
    annotate("text",x=xlim[2]-1,y=quantile(-log10(p),0.9999),label=group1)+
    theme_bw(base_size = 16) +
    geom_vline(xintercept=c(-0.41,0.41),
               lty=4,col="orange",lwd=1)+ # add vetical line
    geom_hline(yintercept = -log10(0.05),
               lty=4,col="orange",lwd=1)+ #add hori line
    labs(x="log2 (Fold change)",
         y="-log10 (p-value)",
         title="Volcano plot")+
    xlim(xlim)+#add changable xlim
    ggrepel::geom_text_repel(
      data = subset(vol, p < p.cutoff&abs(log2(fc))>fc.cutoff),###fc的绝对值大于1
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
      geom_point(aes(color = Significant),size=3) +
      scale_color_manual(values = colour) +
      annotate("text",x=xlim[2]-1,y=quantile(-log10(p),0.9999)-h,label=group2)+
      annotate("text",x=xlim[2]-1,y=quantile(-log10(p),0.9999),label=group1)+
      theme_bw(base_size = 16) +
      geom_vline(xintercept= 0,
                 lty=4,col="orange",lwd=1)+
      geom_hline(yintercept = -log10(0.05),
                 lty=4,col="orange",lwd=1)+ #
      labs(x="log2 (Fold change)",
           y="-log10 (p-value)",
           title="Volcano plot")+
      xlim(xlim)+
      geom_text_repel(
        data = subset(vol, p < p.cutoff),
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
