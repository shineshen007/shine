#' @title CorAnalysis
#' @description a function can do correlation analysis
#' @param number.cex The cex parameter to send to the call to text when
#' writing the correlation coefficients into the plot.
#' @param number.digits indicating the number of decimal digits to be added
#' into the plot. Non-negative integer or NULL, default 4.
#' @param adjust What adjustment for multiple tests should be used
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
CorAnalysis<-function(number.cex = 0.6,number.digits=4,
                      adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")){
  require(corrplot);require(PerformanceAnalytics)
  require(psych)
  data1<-read.csv("data1.csv",check.names = F)
  data2<-read.csv("data2.csv",check.names = F)
  core1<-as.data.frame(data1[,-1])
  core2<-as.data.frame(data2[,-1])
  cor<-corr.test(core1,core2,adjust = adjust)
  #cor plot
  r<-cor[["r"]]
  tiff(file="bicor plot.tiff", width = 1200, height = 1000,res = 56*2)
  col=colorRampPalette(c("navy","white","firebrick3"))
  corrplot(r,tl.col=col(10), tl.srt=45,tl.cex = 0.8,number.cex = number.cex,
           addCoefasPercent = TRUE,cl.lim = c(-1,1),addCoef.col = "black")
  dev.off()

  p<-cor[["p"]]
  tiff(file="pvalue plot.tiff", width = 1200, height = 1000,res = 56*2)
  corrplot(p,tl.col="black", tl.srt=45,tl.cex = 0.8,number.cex = number.cex,
           cl.lim = c(0,1),addCoef.col = "black",number.digits = number.digits)
  dev.off()

  #chart plot
  #tiff(file="chart plot.tiff", width = 1200, height = 1000,res = 56*2)
  #chart.Correlation(core,histogram = TRUE,pch=19)
  #dev.off()
  #heatmap
  #tiff(file="heatmap.tiff", width = 1200, height = 1000,res = 56*2)
  #col = colorRampPalette(c("blue", "white", "red"))(20)
  #heatmap(x = cor, col = col, symm = TRUE)
  #dev.off()
}
