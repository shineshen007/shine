#' @title Scatter3D
#' @description a function to draw 3D scatter plot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
Scatter3D <- function(){
  #library(plot3D)
  data <- read_csv('data.csv')
  data$fc <- log2(data$fc)
  data$p <- -log10(data$p)
  colnames(data)[4] <- 'vip'
  x <- data$fc
  y <- data$p
  z <- data$vip
  scatter3D(x,y,z,bty = 'b2', colkey = FALSE,phi=0,pch=16, cex=1.5,
            ticktype = "detailed",col='DarkMagenta',
            xlab = 'log2(Foldchange)',ylab = '-log10(FDR)',zlab = 'vip')
  text3D(x,y,z,labels=data$compound.name,add = T)
  export::graph2ppt(file='3d.pptx',width=9,height=7)
}

