#' @title OPLSDA
#' @description a function to do OPLSDA
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param group group
#' @param data_position 1
#' @param scale method to scale data
#' @param repeats 200
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
OPLSDA <- function(group='Gout',
                       data_position = 1,
                       scale="ParetoNorm",
                       repeats=200){
  # library(MetaboAnalystR)
  # library(export)
  # library(ggplot2)
  # library(readr)
  # library(magrittr)
  #data pre
  data <- read_csv(dir()[data_position])
  sample.info <- read_csv("sample.info.csv")
  Subject.idx<-sample.info$sample.name[sample.info$class=="Subject"] %>%
    match(.,colnames(data))#get Subject index
  sample <- data[,Subject.idx]
  dm <- cbind(data$name,sample)
  sam_group <- sample.info$group[sample.info$class=="Subject"]
  label <- as.matrix(c('Label',sam_group))
  #colnames(label) <- colnames(dm)
  dma <- cbind(label,t(dm)) %>% t()
  oda <- dma[,-which(dma[1,]==group)]
  #
  fn <- paste('oplsda',unique(unique(sam_group)[- which(unique(sam_group)==group)])[1])
  fn <- paste(fn,unique(unique(sam_group)[- which(unique(sam_group)==group)])[2],sep = '&')
  dir.create(fn)
  setwd(fn)
  write.csv(oda,'oplsda_data.csv',row.names = F)
  #
  path <- getwd()
  data_path <- paste(path,'/oplsda_data.csv',sep = '')
  mSet<-InitDataObjects("pktable", "stat", FALSE)
  mSet<-Read.TextData(mSet,data_path, "colu", "disc");
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet);
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "NULL", scale, ratio=FALSE, ratioNum=20)
  mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
  mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
  mSet<-OPLSR.Anal(mSet, reg=TRUE)
  mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
  mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
  mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", "png", 72, width=NA)
  mSet<-UpdateOPLS.Splot(mSet, "all");
  mSet<-PlotOPLS.Splot(mSet, "opls_splot_1_", "all", "png", 72, width=NA);
  mSet<-UpdateOPLS.Splot(mSet, "custom");
  mSet<-PlotOPLS.Splot(mSet, "opls_splot_2_", "custom", "png", 72, width=NA);
  mSet<-OPLSDA.Permut(mSet, repeats)
  mSet<-PlotOPLS.Permutation(mSet, "opls_perm_1_", "png", 72, width=NA)
  data <- read.csv("oplsda_score.csv")
  colnames(data) <- c('sample','Score','OrthoScore')

  opls <- ggplot(data=data,aes(x= Score,y=OrthoScore,color=mSet$dataSet$cls))+
    geom_point(size=3)+
    scale_color_brewer(palette = 'Set2')+
    #scale_color_discrete(name = "My Legend Title No. 1")+
    stat_ellipse()+
    labs(x=paste("T score [1]", "(", round(100*mSet$analSet$oplsda$modelDF["p1", "R2X"],1), "%)"),
         y=paste("Orthogonal T score [1]", "(", round(100*mSet$analSet$oplsda$modelDF["o1", "R2X"],1), "%)"),
         title="OPLS-DA",
         color = "Group")+
    geom_vline(xintercept=0,
               lty=2,lwd=0.5)+ # add vetical line
    geom_hline(yintercept = 0,
               lty=2,lwd=0.5)+ #add hori line
    theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),#the font size of axis
          axis.title.x = element_text(size = 14),#the font size of axis title
          axis.title.y = element_text(size = 14),#remove ggplot2 background
          legend.title = element_text(size=14),legend.text = element_text(size = 14),
          panel.background = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"))
  FileName <- paste('OPLSDA', levels(mSet[["dataSet"]][["proc.cls"]])[1], sep = " ")
  FileName <- paste(FileName, levels(mSet[["dataSet"]][["proc.cls"]])[2], sep = " ")
  FileName <- paste(FileName, scale, sep = " ")
  FileName <- paste(FileName, ".pptx", sep = "")
  save(opls,file='opls.Rda')
  export::graph2ppt(x=opls,width=9,height=7,file=FileName)
  #permutation
  perm.res <- mSet$analSet$oplsda$perm.res;

  r.vec <- perm.res$r.vec;
  q.vec <- perm.res$q.vec;
  perm.num <- perm.res$perm.num;
  YY<-as.numeric(mSet[["dataSet"]][["proc.cls"]])
  cor <- NULL
  for (i in 1:repeats) {
    Yt <- sample(YY)
    cor[i] <- abs(cor(Yt, YY))
    cat(i); cat(" ")
  }
  qr <- as.data.frame(cbind(r.vec[-1],q.vec[-1],cor))
  #
  lm.r2 <- lm(c(r.vec[1],r.vec[-1])~c(1, cor))
  lm.q2 <- lm(c(q.vec[1],q.vec[-1])~c(1, cor))
  intercept.q2 <- lm.q2$coefficients[1]
  intercept.r2 <- lm.r2$coefficients[1]
  slope.q2 <- lm.q2$coefficients[2]
  slope.r2 <- lm.r2$coefficients[2]
  #
  per <- ggplot2::ggplot(qr)+
    geom_point(aes(x=cor,y=q.vec[-1],color='Q2'),size=3) +
    geom_point(aes(x = cor, y = r.vec[-1],color='R2'),size=3)+
    geom_point(aes(x = 1, y = r.vec[1],color='R2'),size=3)+
    geom_point(aes(x=1,y=q.vec[1],color='Q2'),size=3) +
    scale_color_manual(labels = c(paste("Q2",round(q.vec[1],2), sep = ": "),
                                  paste("R2",round(r.vec[1],2),sep=": ")),
                       values = c("royalblue", "palegreen"),name = "Intercepts") +
    theme_bw(base_size = 16) +
    theme(legend.position = 'bottom',legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),#the font size of axis
          axis.title.x = element_text(size = 14),#the font size of axis title
          axis.title.y = element_text(size = 14))+
    geom_segment(x = 0, y = intercept.q2, xend = 1, yend = q.vec[1], lty = 2)+
    geom_segment(x = 0, y = intercept.r2, xend = 1, yend = r.vec[1], lty =2)+
    labs(x="Correlation",
         y="Values (Q2, R2)",
         title="Permutation test")+
    xlim(0:1)+
    ylim(c(min(c(q.vec[-1],r.vec[-1],q.vec[1],r.vec[1])),1.2*max(c(q.vec[-1],r.vec[-1],q.vec[1],r.vec[1]))))
  save(per,file='permutation.Rda')
  export::graph2ppt(x=per,file=FileName,width=9,height=7,append=T)
  #splot
  sp <- read.csv("oplsda_splot.csv")
  index<-c(1:nrow(sp))
  datas <- cbind(index,sp)
  # Significant_s<- as.factor(ifelse(abs(log2(datas$fc)) > 0.41,
  #                                  ifelse(log2(datas$fc) < -0.41,
  #                                         "Down","Up"),"Not Sig"))
  splot <- ggplot2::ggplot(datas, aes(x = p.1., y = p.corr..1.))+
    geom_point()+
    #scale_color_gradient(low="blue", high="red")+
    #scale_color_brewer(palette = 'Set2')+
    labs(title="Splot of OPLS-DA")+
    geom_vline(xintercept=0,col="orange",
               lty=2,lwd=0.5)+ # add vetical line
    geom_hline(yintercept = 0,col="orange",
               lty=2,lwd=0.5)+ #add hori line
    labs(x='p[1]',
         y='p(corr)[1]',
         title="Splot of OPLS-DA",
         color = "p(corr)[1]")+
    theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),#remove ggplot2 background
          legend.position = 'none',axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),#legend.text = element_text(size = 14),
          axis.title.x = element_text(size = 14),#the font size of axis title
          axis.title.y = element_text(size = 14),
          panel.background = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"))+
    ggrepel::geom_text_repel(
      data = subset(datas, abs(p.corr..1.) > 0.4&abs(p.1.)>max(abs(sp$p.1.))*0.5),
      max.iter = 100000,
      aes(label = X),
      size = 4,
      box.padding = 0.25,
      point.padding = 0.3
    )
  save(splot,file='splot.Rda')
  export::graph2ppt(x=splot,file=FileName,width=9,height=7,append=T)
  setwd('..//')
}
