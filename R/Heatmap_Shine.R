#' @title Heatmap_Shine
#' @description a function can generate heatmap
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param colour colour selection.
#' @param a lowerlimt of colour.
#' @param b toplimit of colour
#' @param c (abs(a)+b)/d
#' @param d default 100
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
Heatmap_Shine <- function(colour = c("green","white","red"),a=-2,b=2,c=0.04,
                          d=100){
  require(pheatmap);require(data.table)
  cat("Import data...\n")
  data <- fread("data.csv")
  data <- setDF(data)
  sample.info <- read.csv("sample.info.csv")
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  sample<-data[,match(sample.name,colnames(data))]

  cat("Draw Heatmap...\n")
  #### heatmap
  x<-sample
  y<-as.matrix(x)
  scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
  }
  y<-scale_rows(y)
  y[y>b]=b
  y[a>y]=a
  bk = unique(c(seq(a,b,c)))
  anno<-data.frame(sample.info[,-c(2:4)],row.names = T)
  tiff(file="heatmap.tiff", width = 1600, height = 1200,res = 56*2)
  hm <- pheatmap::pheatmap(y,color=colorRampPalette(colour)(d),
                           border_color = "grey60",
                           scale = "none",
                           breaks = bk,
                           fontsize=10,
                           fontsize_row=8,
                           fontsize_col=6,
                           annotation=anno)

  dev.off()
}
