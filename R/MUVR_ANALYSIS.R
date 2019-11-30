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
MUVR_ANALYSIS <- function(group = "N",#the group you want to remove
                          data_position = 1,
                          info_position = 2
){

  data <- readr::read_csv(dir()[data_position])
  info <- readr::read_csv(dir()[info_position])
  sample.name<-info$sample.name[info$class=="Subject"]

  sample<-data[,match(sample.name,colnames(data))]%>%
    t(.)
  colnames(sample) <- data$name

  for (i in 1:nrow(sample)) {
    rownames(sample)[i]=info$group[info$sample.name==rownames(sample)[i]]
  }

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
  vmid <- getVIP(classModel, model = 'mid')
  rdmid <- data[match(intersect(row.names(vmid),data$name),data$name),] %>%
    write.csv(., 'MUVR mid.csv', row.names = F)
  vmax <- getVIP(classModel, model = 'max')
  rdmax <- data[match(intersect(row.names(vmax),data$name),data$name),] %>%
    write.csv(., 'MUVR max.csv', row.names = F)
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
