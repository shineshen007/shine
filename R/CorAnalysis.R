#' @title CorAnalysis_Shine
#' @description a function can do correlation analysis
#' @param number.cex The cex parameter to send to the call to text when
#' writing the correlation coefficients into the plot.
#' @param tl.cex font size of label
#' @param number.digits indicating the number of decimal digits to be added
#' into the plot. Non-negative integer or NULL, default 4.
#' @param adjust What adjustment for multiple tests should be used
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @import corrplot
CorAnalysis_Shine <- function(number.cex = 0.6,#number size in cicle
                      tl.cex = 0.6,#font size of label
                      number.digits=4,
                      adjust = c("holm", "bonferroni", "BH", "fdr", "none")){

  data1<-read.csv("data1.csv",check.names = F)
  data2<-read.csv("data2.csv",check.names = F)
  core1<-as.data.frame(data1[,-1])
  core2<-as.data.frame(data2[,-1])
  cor<-psych::corr.test(core1,core2,adjust = adjust)
  #cor plot
  r<-cor[["r"]]
  #unique cor
  pjc <- NULL
  for (i in 1:ncol(r)) {
    padc <- rep(rownames(r)[i],ncol(r))%>%
      #cbind(.,colnames(pa))%>%
      cbind(.,r[i,])
    pjc <- rbind(pjc,padc)
  }
  #png(file="bicor plot.png", width = 1200, height = 1000,res = 56*2)
  col=colorRampPalette(c("navy","white","firebrick3"))
  cp <- corrplot::corrplot(r,tl.col="black", tl.srt=45,tl.cex = tl.cex,number.cex = number.cex,
                     addCoefasPercent = TRUE,cl.lim = c(-1,1),addCoef.col = "black",col = col(10))
  save(cp,file = 'correlation.Rda')
  export::graph2ppt(x=cp,file='correlation.pptx',height=7,width=9)
  #dev.off()

  p<-cor[["p"]]
  #png(file="pvalue plot.png", width = 1200, height = 1000,res = 56*2)
  cpp <- corrplot::corrplot(p,tl.col="black", tl.srt=45,tl.cex = tl.cex,number.cex = number.cex,
                     cl.lim = c(0,1),addCoef.col = "black",number.digits = number.digits)
  save(cpp,file = 'correlation_p.Rda')
  export::graph2ppt(x=cpp,file='correlation.pptx',height=7,width=9,append = TRUE)
  #dev.off()
  #unique p
  pj <- NULL
  for (i in 1:ncol(p)) {
    pad <- rep(rownames(p)[i],ncol(p))%>%
      #cbind(.,colnames(pa))%>%
      cbind(.,p[i,])
    pj <- rbind(pj,pad)
  }

  #
  result <- cbind(pjc,pj)#%>%
  result <-  result[, -3]

  colnames(result) = c('M2',"cor",'pval')
  write.csv(result,'result.csv')
  #rm p>0.05&pcor=1 generate data for cyto
  rst <- read.csv('result.csv')
  rs <- rst[-which(rst$pval>0.05),]
  rs <- rs[-which(rs$cor==1),]
  colnames(rs)[1] <- 'M1'
  dup <- which(duplicated(rs$cor))
  rsu <- rs[-dup,]
  xlsx::write.xlsx(rsu,'cor network.xlsx',row.names = F)
  #chart plot
  #png(file="chart plot.png", width = 1200, height = 1000,res = 56*2)
  #chart.Correlation(core,histogram = TRUE,pch=19)
  #dev.off()
  #heatmap
  #png(file="heatmap.png", width = 1200, height = 1000,res = 56*2)
  #col = colorRampPalette(c("blue", "white", "red"))(20)
  #heatmap(x = cor, col = col, symm = TRUE)
  #dev.off()
}
