#' @title ArrangeData
#' @description a function to arrange data for svr
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
ArrangeData <- function(){
  data <- data.table::fread("peak.table.csv")
  data<- data.table::setDF(data)

  data<-data[,-c(3:4,6:10)]#remove redundancy columns
  colnames(data)[2] <- 'mz'
  colnames(data)[3] <- 'rt'
  write.csv(data,"data for svr.csv",row.names = F)

}
