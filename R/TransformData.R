#' @title TransformData
#' @description a function to TransformData
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data_name data_position
#' @param name compound.name
#' @param boxplot_data default is TRUE
#' @return  All the results can be got form other functions and instruction.
#' @export
TransformData <- function(
  data_name = "data.csv",
  name='compound.name',
  boxplot_data = TRUE
){
  data <- readr::read_csv(data_name) %>% as.data.frame()
  info <- readr::read_csv('sample.info.csv')
  sample.name<-info$sample.name[info$class=="Subject"]

  sample<-data[,match(intersect(colnames(data),sample.name),colnames(data)),]%>%
    t(.)
  colnames(sample) <- data[,name]
  rownames(sample)=info$group[match(rownames(sample),info$sample.name)]
  if(boxplot_data){
    #data for boxplot
    df <- NULL
    for (i in 1:ncol(sample)) {
      ads <- cbind(rownames(sample),scale(sample[,i]),rep(colnames(sample)[i],nrow(sample)))
      adf <- NULL
      ds <- rbind(adf,ads)
      df <- rbind(df,ds)
    }
    colnames(df) <- c('group','abundance','metabolites')
    write.csv(df,'data for boxplot.csv',row.names = F)
  }

  write.csv(sample,'data for roc.csv')

}



