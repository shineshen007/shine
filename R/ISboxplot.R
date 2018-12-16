#' @title ISboxplot
#' @description a function can generate internal standard boxplot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param sample a numeric object show the number of sample.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
ISboxplot <- function(sample = NULL){

  data <- read.csv("IS data.csv",check.names = F)
  mz <- c(rep(c(273.0762,271.0616),sample))
  data1 <- stack(data)
  data1 <- cbind(data1,mz)
  tiff(filename = "boxplot of IS.tiff", width = 800, height = 600,res = 56*2)
  gplot <- ggpubr::ggboxplot(data = data1,x="mz",y="scale(values)",color = "mz",
                  palette = "jco",add = "jitter",xlab = "Internal Standard",
                  ylab = "Peak Area",title = "Boxplot of Internal Standard")
  plot(gplot)
  dev.off()
}
