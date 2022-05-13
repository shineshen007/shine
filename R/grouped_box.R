#' @title Grouped_BoxVilion
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param levels order of group
#' @param palette palette
#' @param stat_type parametric
#' @param plot.type box
#' @param data_name data.csv
#' @param zscore default TRUE
#' @param pairwise.display significant
#' @param p.adjust.method fdr
#' @param plot_nrow plot_nrow
#' @param info_name info_name
#' @param ylab name
#' @return  All the results can be got form other functions and instruction.
#' @export
Grouped_BoxVilion <- function(
  stat_type = 'nonparametric',
  data_name = 'data.csv',
  info_name = 'info.csv',
  plot.type = "box",
  zscore=TRUE,
  palette = "nrc_npg",#ggsci
  plot_nrow=2,
  pairwise.display = "significant",
  p.adjust.method = "fdr",
  ylab = "Relative Abundance (log2)",
  levels =c ("M1","M2","M3")
){

  pacman::p_load(PMCMRplus,ggsci,ggstatsplot)
  #if(!file.exists("data for boxplot.csv")){
  #data_name = "3mt_rna.csv"
  data <- read_csv(data_name)
  colnames(data)[1]='name'
  info <- read_csv(info_name)
  sample.name<-info$sample.name[info$class=="Subject"]

  sample<-data[,match(sample.name,colnames(data))]%>%
    t(.)
  colnames(sample) <- data$name

  if(zscore){
    #data for boxplot
    for (i in 1:nrow(sample)) {
      rownames(sample)[i]=info$group[which(info$sample.name==rownames(sample)[i])]
    }

    nc <- ncol(sample)
    df <- NULL
    #log2
    for (i in 1:nc) {
      ads <- cbind(rownames(sample),scale(sample[,i]),rep(colnames(sample)[i],nrow(sample)))
      adf <- NULL
      ds <- rbind(adf,ads)
      df <- rbind(df,ds)
    }
    colnames(df) <- c('group','abundance','metabolites')
    write.csv(df,'data for boxplot.csv',row.names = F)
  }else{
    #data for boxplot
    for (i in 1:nrow(sample)) {
      rownames(sample)[i]=info$group[which(info$sample.name==rownames(sample)[i])]
    }

    nc <- ncol(sample)
    df <- NULL
    #log2
    for (i in 1:nc) {
      ads <- cbind(rownames(sample),sample[,i],rep(colnames(sample)[i],nrow(sample)))
      adf <- NULL
      ds <- rbind(adf,ads)
      df <- rbind(df,ds)
    }
    colnames(df) <- c('group','abundance','metabolites')
    write.csv(df,'data for boxplot.csv',row.names = F)
  }
  #}
  #rm(list = ls())

  data <- read.csv("data for boxplot.csv")
  data$group <- factor(data$group, levels = levels)
  data$metabolites <- fct_inorder(data$metabolites)


  groupbox <- grouped_ggbetweenstats(data = data,x=group,y = abundance,
                                     results.subtitle = F,
                                     grouping.var = metabolites,
                                     #pairwise.comparisons = T,
                                     #title = 'enyname of 3MT',
                                     plot.type = plot.type,
                                     type = stat_type,
                                     plotgrid.args = list(nrow=plot_nrow),
                                     ylab = ylab,
                                     pairwise.display = pairwise.display, ## display only significant pairwise comparisons
                                     p.adjust.method = p.adjust.method,
                                     package = "ggsci",
                                     palette = palette)
  # theme(axis.text.x=element_text(size=font_size),
  #       legend.title=element_blank(),
  #       legend.text = element_text(size = font_size),
  #       axis.text.y = element_text(size = font_size),
  #       axis.title.x = element_text(size = 0),#the font size of axis title
  #       axis.title.y = element_text(size = font_size))
  save(groupbox,file = 'group_box.Rda')
  ggsave('group_box.pdf')

}
