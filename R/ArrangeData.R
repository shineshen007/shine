#' @title ArrangeData
#' @title ArrangeData
#' @description a function to arrange data for svr
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
ArrangeData <- function(){
  require(data.table)
  data <- fread("peaktable.csv")
  data<- setDF(data)
  info <- fread("sample.info.csv")
  info<-setDF(info)
  info<-info[,-1]

  data<-data[,-c(2:4,6,7,9:13)]#remove redundancy column
  colnames(data)[2] <- 'mz'
  colnames(data)[3] <- 'rt'
  write.csv(data,"data for svr.csv",row.names = F)
  colnames<-as.data.frame(colnames(data))
  colnames<-as.data.frame(colnames[-c(1:3),])
  colnames(colnames)[1] <- 'sample.name'
  write.csv(colnames,"sample.name.csv",row.names = F)
  info<-cbind(colnames,info)
  write.csv(info,"info.csv",row.names = F)

}
