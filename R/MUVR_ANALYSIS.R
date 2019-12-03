#' @title MUVR_ANALYSIS
#' @description a function to do MUVR_ANALYSIS
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param group group
#' @param data_position 1
#' @param info_position 2
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
MUVR_ANALYSIS <- function(group = "Normal",#the group you want to remove
                          data_position = 1,
                          info_position = 2
){

  data <- readr::read_csv(dir()[data_position])
  info <- readr::read_csv(dir()[info_position])
  sample.name<-info$sample.name[info$class=="Subject"]

  sample<-data[,match(intersect(sample.name,colnames(data)),colnames(data))]%>%
    t(.)
  colnames(sample) <- data$name
  rownames(sample)=info$group[match(rownames(sample),info$sample.name)]
  #data for roc
  s <- which(rownames(sample) == group)
  as <- sample[-s,]
  write.csv(as,'ls.csv')
  roc <- read_csv('ls.csv')
  colnames(roc)[1] <- 'group'
  unlink('ls.csv')
  #MUVR
  # library(MUVR)
  # library(parallel)
  # library(doParallel)
  data_MUVR <- as.data.frame(roc)
  # Check number of observations per class
  X <- data_MUVR[, -1]
  Y <- data_MUVR[,1]
  # Set method parameters
  nCore = parallel::detectCores() - 1   # Number of processor threads to use
  nRep = nCore              # Number of MUVR repetitions
  nOuter = 8                # Number of outer cross-validation segments
  varRatio = 0.8            # Proportion of variables kept per iteration
  method = 'RF'             # Selected core modelling algorithm
  # Set up parallel processing using doParallel

  cl = parallel::makeCluster(nCore)
  registerDoParallel(cl)
  # Perform modelling
  set.seed(1234)
  classModel = MUVR(
    X = X,
    Y = Y,
    nRep = nRep,
    nOuter = nOuter,
    varRatio = varRatio,
    method = method
  )
  save(classModel,file = 'MUVR.RData')
    # Stop parallel processing
  stopCluster(cl)
  # Examine model performance and output
  # classModel$miss                   # Number of misclassifications for min, mid and max models
  #
  # classModel$nVar                   # Number of variables for min, mid and max models

  a <-cbind(Y, classModel$yClass)    # Actual class side-by-side with min, mid and max predictions
  plotVAL(classModel)

  export::graph2ppt(file = 'MUVR.pptx',
                    height = 7,
                    width = 9)
  v <- getVIP(classModel, model = 'min')
  rdmin <- data[match(intersect(row.names(v),data$name),data$name),] %>%
    write.csv(., 'MUVR min.csv', row.names = F)
  vmid <- getVIP(classModel, model = 'mid')
  rdmid <- data[match(intersect(row.names(vmid),data$name),data$name),] %>%
    write.csv(., 'MUVR mid.csv', row.names = F)
  vmax <- getVIP(classModel, model = 'max')
  rdmax <- data[match(intersect(row.names(vmax),data$name),data$name),] %>%
    write.csv(., 'MUVR max.csv', row.names = F)
  #
  #
  iplot <- ggplot2::ggplot(v, aes(x = rank , y = reorder(name,rank)))+
    geom_point(aes(color = name),size=3) +
    labs(title="metabolites importance plot")+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),#the font size of axis title
          axis.title.y = element_text(size = 14))+
    xlab('importance')+
    theme(legend.position="none")
  save(iplot,file = 'importance.Rda')
  #ggsave('metabolites importance plot.png', width = 12, height = 8)
  export::graph2ppt(x=iplot,file='MUVR.pptx',height=7,width=9,append=TRUE)
  #
  ia <- NULL
  for (i in 1:nrow(v)) {
    r <- grep(row.names(v)[i], colnames(data_MUVR))
    ia <- c(ia, r)
  }
  dr <- data_MUVR[, c(1, ia)]
  colnames(dr)[1] <- 'group'
  write.csv(dr, 'data for roc.csv', row.names = F)
  #get roc data
  rd <- data[match(intersect(colnames(dr)[-1],data$name),data$name),]
  xlsx::write.xlsx(rd,'roc data.xlsx')
}
