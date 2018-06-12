#' @title Heatmap
#' @description a function can generate heatmap
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns.
#' @param sample.info a dataframe include sample.name,injection.order,
#' class,batch and group columns.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
Heatmap <- function(data = NULL,sample.info = NULL){
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
  bk = unique(c(seq(-2,2, 0.04)))
  anno<-data.frame(sample.info[,-c(2:4)],row.names = T)
  png(file="heatmap.png", width = 1600, height = 1200,res = 56*2)
  hm <- pheatmap::pheatmap(y,color=colorRampPalette(c("green","white","red"))(100),
                           border_color=NA,
                           scale = "row",
                           breaks = bk,
                           fontsize=10,
                           fontsize_row=8,
                           fontsize_col=6,
                           annotation=anno)

  dev.off()
}
