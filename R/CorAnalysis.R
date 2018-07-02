#' @title CorAnalysis
#' @description a function can do correlation analysis
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data a dataframe to do correlation analysis.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
CorAnalysis<-function(data){
  require(corrplot);require(PerformanceAnalytics)
  data<-read.csv("data.csv")
  core<-data[,-1]
  cor<-cor(core)
  #co plot
  png(file="cor plot.png", width = 1200, height = 1000,res = 56*2)
  corrplot(cor,type="upper", order="hclust", tl.col="black", tl.srt=45,tl.cex = 0.8)
  dev.off()
  #chart plot
  png(file="chart plot.png", width = 1200, height = 1000,res = 56*2)
  chart.Correlation(core,histogram = TRUE,pch=19)
  dev.off()
  #heatmap
  png(file="heatmap.png", width = 1200, height = 1000,res = 56*2)
  col = colorRampPalette(c("blue", "white", "red"))(20)
  heatmap(x = cor, col = col, symm = TRUE)
  dev.off()
}
