#' @title CorAnalysis
#' @description a function can do correlation analysis
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
CorAnalysis<-function(){
  require(corrplot);require(PerformanceAnalytics)
  data1<-read.csv("data1.csv")
  data2<-read.csv("data2.csv")
  core1<-data1[,-1]
  core2<-data2[,-1]
  cor<-cor(core1,core2)
  #cor plot
  png(file="bicor plot.png", width = 1200, height = 1000,res = 56*2)
  corrplot(cor,tl.col="black", tl.srt=45,tl.cex = 0.8,
           addCoefasPercent = TRUE,cl.lim = c(-100,100),addCoef.col = "black")
  dev.off()
  cor1<-cor(core1)
  png(file="cor plot.png", width = 1200, height = 1000,res = 56*2)
  corrplot(cor1,type="upper", order="hclust", tl.col="black", tl.srt=45,tl.cex = 0.8,
           addCoefasPercent = TRUE,cl.lim = c(-100,100),addCoef.col = "black")
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
