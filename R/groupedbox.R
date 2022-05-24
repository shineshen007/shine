#' @title Grouped_Box
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param levels order of group
#' @param font_size font_size
#' @param zscore default TRUE
#' @param title name
#' @param ylab name
#' @param xlab_angle 45
#' @param data_name data_name
#' @param info_name info_name
#' @param variable_base variable_base
#' @param group_base group_base
#' @return  All the results can be got form other functions and instruction.
#' @export
Grouped_Box <- function(
  font_size = 15,
  levels = c("M1","M2","M3"),
  zscore=TRUE,
  xlab_angle = 45,
  title = NULL,
  ylab = "Relative Abundance (log2)",
  info_name = 'info.csv',
  data_name = "data.csv",
  group_base = TRUE,
  variable_base = FALSE
){
  pacman::p_load(ggsci,ggpubr)
  #if(!file.exists("data for boxplot.csv")){

  data <- read_csv(data_name)
  colnames(data)[1]='name'
  info <- read_csv(info_name)
  sample.name<-info$sample.name[info$class=="Subject"]

  sample<-data[,match(intersect(sample.name,colnames(data)),colnames(data))]%>%
    t(.)
  colnames(sample) <- data$name

  if(zscore){
    #data for boxplot
    for (i in 1:nrow(sample)) {
      rownames(sample)[i]=info$group[which(info$sample.name==rownames(sample)[i])]
    }

    nc <- ncol(sample)
    df <- NULL
    #
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

  data <- read.csv("data for boxplot.csv")
  data$group <- factor(data$group, levels = levels)
  data$metabolites <- fct_inorder(data$metabolites)
  my_comparisons <- list(c(levels[1],levels[2]), c(levels[2], levels[3]),
                         c(levels[1], levels[3]))
  if(variable_base){
    s <- ggplot(data,aes(x=metabolites,y = abundance))+
      geom_boxplot(aes(fill = group))+
      theme_bw()+
      scale_fill_npg()+
      labs(title=title)+
      ylab(ylab)+
      theme(axis.text.x=element_text(size=font_size),#angle=xlab_angle,hjust = 1,vjust = 1,
            plot.title = element_text(size=font_size,hjust = 0.5),
            legend.title=element_blank(),
            legend.position = 'none',
            legend.text = element_text(size = font_size),
            axis.text.y = element_text(size = font_size),
            axis.title.x = element_text(size = 0),#the font size of axis title
            axis.title.y = element_text(size = 0)
      )+
      stat_compare_means(aes(group = group),label = "p.signif", #comparisons = my_comparisons,
                         label.x = 1.5)
    save(s,file = 'variable_base_Box.rda')
    ggsave('variable_base_Box.pdf',width = 9,height = 6)
  }
  if(group_base){
    s <- ggplot(data,aes(x=group,y = abundance))+
      geom_boxplot(aes(fill = group))+
      theme_bw()+
      scale_fill_npg()+
      labs(title=title)+
      ylab(ylab)+
      theme(axis.text.x=element_text(size=font_size),
            plot.title = element_text(size=font_size,hjust = 0.5),
            legend.title=element_blank(),
            legend.position = 'none',
            legend.text = element_text(size = font_size),
            axis.text.y = element_text(size = font_size),
            axis.title.x = element_text(size = 0),#the font size of axis title
            axis.title.y = element_text(size = 0)
      )+
      stat_compare_means(comparisons = my_comparisons)+ #添加成对p值
      stat_compare_means(label.y = max(data$abundance)+5)
    save(s,file = 'group_base_Box.rda')
    ggsave('group_base_Box.pdf',width = 9,height = 6)
  }
  #save(s,file = 'Grouped_Box.Rda')
  #export::graph2ppt(x=s,file='violinplot.pptx',height=12,width=ppt_width)
}




