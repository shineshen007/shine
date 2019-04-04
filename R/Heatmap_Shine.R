#' @title Heatmap_Shine
#' @description a function can generate heatmap
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param colour colour selection.
#' @param a lowerlimt of colour.
#' @param b toplimit of colour
#' @param c (abs(a)+b)/d
#' @param d default 100
#' @param group group info.
#' @param scale_row scale or not
#' @param size_row the font size of row
#' @param size_col the font size of col
#' @param fontsize the font size
#' @param border show boder or not
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
Heatmap_Shine <- function(colour = c("green","white","red"),
                          a=-2,#lower limit
                          b=2,#upper limit
                          c=0.04,#(a+b)/d
                          d=100,#interval
                          scale_row = TRUE,
                          size_row=10,
                          size_col=8,
                          fontsize=10,
                          border = NA,
                          group = c("case","control")
){

  cat("Import data...\n")
  data <- data.table::fread("data pathway.csv")
  data <- data.table::setDF(data)
  sample.info <- read.csv("sample.info.csv")

  case.name<-sample.info$sample.name[sample.info$group==group[1]] #get case name
  control.name<-sample.info$sample.name[sample.info$group==group[2]]
  qc.name<-sample.info$sample.name[sample.info$class=="QC"]#get qc index
  qc.idx<-match(qc.name,colnames(data))#get qc index
  data <- data[,-qc.idx]

  case<-data[,match(case.name,colnames(data))]
  control<-data[,match(control.name,colnames(data))]
  sample <- cbind(case,control)
  cat("Draw Heatmap...\n")
  #### heatmap
  x<-sample
  y<-as.matrix(x)
  scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
  }
  if(scale_row){
    y<-scale_rows(y)
  }

  y[y>b]=b
  y[a>y]=a
  bk = unique(c(seq(a,b,c)))
  cg <- c(rep(group[1],length(case.name)))
  cog <- c(rep(group[2],length(control.name)))
  cac <- c(cg,cog)
  anno <- as.data.frame(cac,row.names = T)
  rownames(anno) <- colnames(y)
  colnames(anno) <- "group"
  hm <- pheatmap::pheatmap(y,color=colorRampPalette(colour)(d),
                           border_color = border,
                           scale = "none",
                           breaks = bk,
                           fontsize=fontsize,
                           cluster_cols = F,
                           fontsize_row=size_row,
                           fontsize_col=size_col,
                           annotation=anno,
                           filename = "heatmap.png"
  )
  export::graph2ppt(x=hm,file='heatmap.pptx',height=7,width=9)

}
