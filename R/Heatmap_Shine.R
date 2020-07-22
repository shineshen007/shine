#' @title Heatmap
#' @description a function can generate heatmap
#' @author Xia Shen
#' \email{shenxia@shanghaitech.edu.cn}
#' @param colour colour selection.
#' @param a lowerlimt of colour.
#' @param b toplimit of colour
#' @param cluster_cols default is FALSE
#' @param data data name
#' @param group group info.
#' @param scale_row scale or not
#' @param fontsize the font size
#' @param QC QC exist
#' @param width 12
#' @param height 7
#' @param show_colnames show_colnames
#' @param show_rownames show_rownames
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' #the data must be named as 'data.csv' and 'sample.info.csv'
#' Heatmap()
#' }
Heatmap_Complex <- function(data = 'data.csv',
                            colour = c("navy", "white", "firebrick3"),
                            a=-2,#lower limit
                            b=2,#upper limit
                            #c=0.04,(a+b)/d
                            #d=100,#interval
                            width = 12,
                            height =7,
                            QC=TRUE,
                            scale_row = TRUE,
                            #size_row=10,
                            #size_col=8,
                            fontsize=10,
                            #cellwidth=NA,
                            #cellheight=NA,
                            #border = NA,
                            show_colnames=FALSE,
                            show_rownames=FALSE,
                            cluster_cols = FALSE,
                            group = c("case","control")
){
  pacman::p_load(ComplexHeatmap,ggplotify)
  ##save parameters
  heatmap.parameters <- c(
    paste(colour, collapse = ","),
    paste(group, collapse = ","),
    a,
    b,
    #d,
    scale_row,
    #size_col,
    fontsize,
    #cellwidth,
    #cellheight,
    #border,
    show_colnames,
    show_rownames,
    cluster_cols
  )
  heatmap.parameters <- data.frame(c(
    'colour',
    'group',
    'a',
    'b',
    #'d',
    'scale_row',
    #'size_col',
    'fontsize',
    #'cellwidth',
    #'cellheight',
    #'border',
    'show_colnames',
    'show_rownames',
    'cluster_cols'
  ),
  heatmap.parameters, stringsAsFactors = FALSE)

  heatmap <- rbind(heatmap.parameters,c("Version", "ezMet"))
  colnames(heatmap) <- c('parameter', 'value')
  write.csv(heatmap,"heatmap.parameters.csv",row.names = F)
  #
  pacman::p_load(data.table,ComplexHeatmap,circlize)
  #
  cat("Import data...\n")
  data <- data.table::fread(data)
  data <- data.table::setDF(data)
  sample.info <- read.csv("sample.info.csv")

  case.name<-sample.info$sample.name[sample.info$group==group[1]] #get case name
  control.name<-sample.info$sample.name[sample.info$group==group[2]]
  if(QC){
    qc.name<-sample.info$sample.name[sample.info$class=="QC"]#get qc index
    qc.idx<-match(qc.name,colnames(data))#get qc index
    data <- data[,-qc.idx]
  }
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
  # cc <- (abs(a)+b)/d
  # bk = unique(c(seq(a,b,cc)))
  cg <- c(rep(group[1],length(case.name)))
  cog <- c(rep(group[2],length(control.name)))
  cac <- c(cg,cog)
  anno <- as.data.frame(cac,row.names = T)
  rownames(anno) <- colnames(y)
  colnames(anno) <- "group"
  rownames(y) <- data$compound.name
  #y <- log2(y)
  set.seed(10)
  annotation1 = HeatmapAnnotation(df = anno)
  mycols <- colorRamp2(breaks = c(a, 0, b),
                       colors = colour)
  h1 <- ComplexHeatmap::Heatmap(y,col=mycols,
                                top_annotation = annotation1,
                                row_names_gp = gpar(fontsize = 9) ,# Text size for row names
                                border = NA,cluster_columns = cluster_cols,
                                show_column_names=show_colnames,show_row_names = show_rownames,
                                heatmap_legend_param = list(
                                  title= NULL, #title_position = "topcenter",
                                  legend_height=unit(8,"cm"), legend_direction="vertical")#,
                                #row_title = "lipid", column_title = "sample"
  ) %>% as.ggplot()
  ggsave('heatmap_complex.pdf',width = width,height = height)
  save(h1,file = 'heatmap.Rda')
  if(nrow(y)<100&ncol(y)<100){
    export::graph2ppt(x=h1,file='heatmap.pptx',height=7,width=9)
  }
  #

}
