#' @title Boxplot_Shine
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
Boxplot_Shine <- function(){
  data <- read.csv("data.csv")
  mn <- ncol(data)
  da <- data.frame(data[,1])
  for (i in 2:mn) {
    dao <- scale(data[,i])
    da <- cbind(da,dao)
  }
  colnames(da) <- colnames(data)
  for (i in 2:mn) {
    s <- ggboxplot(data = da,x='group',y = colnames(da)[i],color = "group",
                   palette = "jco",add = "jitter",size = 1 #,ylab = ylab
    )+
      stat_compare_means()
    plot(s)
    export::graph2ppt(file='boxplot.pptx',height=7,width=9,append = TRUE)
  }
}



