#' @title ISboxplot
#' @description a function can generate internal standard boxplot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param sample a numeric object show the number of sample.
#' @param width the width of plot
#' @param height the height of plot
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
ISboxplot <- function(sample = NULL,width=10,height=6){

  data <- read.csv("IS data.csv",check.names = F)
  RSD <- function(){
    qc_rsd <- data
    rsd_qc <- sapply(seq(nrow(qc_rsd)), function(i){
      SD <- sd(qc_rsd[i,],na.rm = TRUE)
      MEAN <- sum(qc_rsd[i,],na.rm = TRUE)/ncol(data)
      rsd_qc <- SD/MEAN
    })
  }
  rsd <- RSD()
  mz <- c(rep(c(273.0762,271.0616),sample))
  data1 <- stack(data)
  data1 <- cbind(data1,mz)

  gplot <- ggpubr::ggboxplot(data = data1,x="mz",y="scale(values)",color = "mz",
                             palette = "jco",add = "jitter",xlab = "Internal Standard",
                             ylab = "Peak Area",title = "Boxplot of Internal Standard")+
    annotate("text", x="273.0762", y=1.5, label=round(rsd[1],3), colour='black', size=4)+
    annotate("text", x="271.0616", y=1.5, label=round(rsd[2],3), colour='black', size=4)+
    ggsave("IS box plot.png",width=width,height=height)

}
