#' @title Permutation_Shine
#' @description a function to do permutation test.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  get data you want.
#' @param repeats the repeats to do,default is 200
#' @param ncomp the components of pls
#' @param group the group of the data
#' @param scale scale method
#' @export
Permutation_Shine <- function(repeats = 200, ncomp = 2,
                              group = c('case','control'),
                              scale = c('UV','pareto')
) {

  cat("Import data...\n")
  data <- data.table::fread("data.csv")
  data <- data.table::setDF(data)
  sample.info <- readr::read_csv("sample.info.csv")
  class<- sample.info[,"group"]

  group1.index <- which(class == group[1])
  group2.index <- which(class == group[2])
  sample.info<-sample.info[c(group1.index,group2.index),]
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  sample<-data[,match(sample.name,colnames(data))]

  sample.index <- which(sample.info$class=="Subject")
  datat<-sample
  #z-score
  if (scale=="pareto") {
    datatm<-apply(datat,2, function(x) {
      (x-mean(x))/sqrt(sd(x))})}
  if (scale=="UV") {
    datatm<-apply(datat,2, function(x) {
      (x-mean(x))/sd(x)})}

  XXt<-t(datatm)
  group_pls<-as.data.frame(sample.info$group)
  YY<-as.numeric(group_pls[sample.index,])
  pls <- plsdepot::plsreg1(XXt, YY, comps = ncomp)

  Q2 <- pls$Q2[ncomp,5]
  R2 <- sum(pls$R2)

  ##begin repeat
  q2 <- NULL
  r2 <- NULL
  cor <- NULL
  for (i in 1:repeats) {
    Yt <- sample(YY)
    temp.pls <- plsdepot::plsreg1(XXt,Yt,comps = ncomp)
    q2[i] <- temp.pls$Q2[ncomp,5]
    r2[i] <- sum(temp.pls$R2)
    cor[i] <- abs(cor(Yt, YY))
    cat(i); cat(" ")
  }
  qr <- as.data.frame(cbind(q2,r2,cor))
  #
  lm.r2 <- lm(c(R2,r2)~c(1, cor))
  lm.q2 <- lm(c(Q2,q2)~c(1, cor))
  intercept.q2 <- lm.q2$coefficients[1]
  intercept.r2 <- lm.r2$coefficients[1]
  slope.q2 <- lm.q2$coefficients[2]
  slope.r2 <- lm.r2$coefficients[2]
  #
  per <- ggplot2::ggplot(qr)+
    geom_point(aes(x=cor,y=q2,color='Q2'),size=3) +
    geom_point(aes(x = cor, y = r2,color='R2'),size=3)+
    geom_point(aes(x = 1, y = R2,color='R2'),size=3)+
    geom_point(aes(x=1,y=Q2,color='Q2'),size=3) +
    scale_color_manual(labels = c(paste("Q2",round(intercept.q2,2), sep = ": "), paste("R2",round(intercept.r2,2),sep=": ")),
                       values = c("royalblue", "palegreen"),name = "Intercepts") +
    theme_bw(base_size = 16) +
    theme(legend.position = 'bottom',legend.text = element_text(size = 16),
          axis.text.x=element_text(size = 16),
          axis.text.y=element_text(size = 16))+
    geom_segment(x = 0, y = intercept.q2, xend = 1, yend = Q2, lty = 2)+
    geom_segment(x = 0, y = intercept.r2, xend = 1, yend = R2, lty =2)+
    labs(x="Correlation",
         y="Values (Q2, R2)",
         title="Permutation test")+
    xlim(0:1)+
    ylim(c(min(c(q2,r2,Q2,R2)),1.2*max(c(q2,r2,Q2,R2))))
  FileName <- paste('permutation test', group[1], sep = " ")
  FileName <- paste(FileName, group[2], sep = "")
  FileName <- paste(FileName, scale, sep = " ")
  FileName <- paste(FileName, ".pptx", sep = "")
  export::graph2ppt(x=per,file=FileName,width=9,height=7)
}


