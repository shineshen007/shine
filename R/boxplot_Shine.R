#' @title Boxplot_Shine
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param levels order of group
#' @param font_size font_size
#' @param stat_type parametric
#' @param plot.type box
#' @param data_name data.csv
#' @param pairwise.display significant
#' @param p.adjust.method fdr
#' @param ylab name
#' @return  All the results can be got form other functions and instruction.
#' @export
Boxplot_Shine <- function(font_size = 20,
                          stat_type = 'parametric',
                          data_name = 'data.csv',
                          plot.type = "box",
                          pairwise.display = "significant",
                          p.adjust.method = "fdr",
                          ylab = "log2 metabolites abundance",
                          levels =c ("M1","M2","M3")

){
  pacman::p_load(PMCMRplus,ggsci,ggstatsplot,readr)
  data <- read_csv(data_name)
  info <- read_csv(dir()[grep('info.csv',dir())][1])#读取路径中info.csv的文件，如有多个读取第一个
  if(!file.exists('box')){
    dir.create('box')
  }
  file.copy(dir()[grep('info.csv',dir())][1],"box",overwrite = T)
  setwd('box')
  if(!file.exists("data for boxplot.csv")){

    sample.name<-info$sample.name[info$class=="Subject"]
    sample<-data[,match(sample.name,colnames(data))]%>%
      t(.)
    if(nrow(sample)<1|ncol(sample)<1){
      stop("different sample.name in data and info.\n")
    }
    colnames(data)[1]='name'
    colnames(sample) <- data$name
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
  ###
  #rm(list = ls())
  data <- read.csv("data for boxplot.csv")
  #data <- data[1:654,]
  me=unique(data$metabolites)
  for (i in 1:length(me)) {
    cat(i);cat(" ")
    dd <- data[which(data$metabolites==me[i]),]
    dd$group <- factor(dd$group, levels = levels)
    if(stat_type=='nonparametric'){
      pv <- kwAllPairsDunnTest(abundance ~ group, data = dd,p.adjust.method = p.adjust.method)
    }
    if(stat_type=='parametric'){
      pv <- gamesHowellTest(abundance ~ group, data = dd,p.adjust.method = p.adjust.method)
    }

    if(!is.na(any(pv[["p.value"]]<0.05))){
      if(any(pv[["p.value"]]<0.05)){
        #dd$metabolites <- fct_inorder(dd$metabolites)
        fn <- paste0(me[i],'.pdf')
        s <- ggbetweenstats(data = dd,x='group',y = 'abundance',#color = "group",
                            results.subtitle = F,
                            title = me[i],
                            pairwise.comparisons = TRUE,
                            plot.type = plot.type,
                            type = stat_type,
                            #xlab = "Continent",
                            ylab = ylab,
                            pairwise.display = pairwise.display, ## display only significant pairwise comparisons
                            p.adjust.method = p.adjust.method,
                            package = "ggsci",
                            palette = "nrc_npg")+
          theme(axis.text.x=element_text(size=font_size),
                legend.title=element_blank(),
                legend.text = element_text(size = font_size),
                axis.text.y = element_text(size = font_size),
                axis.title.x = element_text(size = 0),#the font size of axis title
                axis.title.y = element_text(size = font_size)
          )
        ggsave(fn)
      }
    }
  }
  #data <- read_csv("..//enrich_score.csv")
  #ss <- data[match(gsub('.pdf','',dir()[grep('.pdf',dir())]),data$...1),]
  setwd('..//')
  #write_csv(ss,'gsva_score.csv')
}




