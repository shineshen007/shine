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
  data1<-read.csv("data1.csv")
  data2<-read.csv("data2.csv")
  core1<-data1[,-1]
  core2<-data2[,-1]
  cor<-corr.test(core1,core2,adjust = adjust)
  #cor plot
  r<-cor[["r"]]
  png(file="bicor plot.png", width = 1200, height = 1000,res = 56*2)
  corrplot(r,tl.col="black", tl.srt=45,tl.cex = 0.8,number.cex = number.cex,
           addCoefasPercent = TRUE,cl.lim = c(-100,100),addCoef.col = "black")
  dev.off()

  p<-cor[["p"]]
  png(file="pvalue plot.png", width = 1200, height = 1000,res = 56*2)
  corrplot(p,tl.col="black", tl.srt=45,tl.cex = 0.8,number.cex = number.cex,
           cl.lim = c(0,1),addCoef.col = "black",number.digits = number.digits)
  dev.off()

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
