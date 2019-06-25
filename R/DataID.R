#' @title DataID
#' @description a function to arrange data
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
DataID <- function(){

  data <- data.table::fread("peak.table.after.data.cleaning.csv")
  info <- data.table::fread("info.csv")
  info<-data.table::setDF(info)
  data<- data.table::setDF(data)
  data<-data[,-1]
  info<-info[,-c(2:4)]
  write.csv(data,"data for dna.csv",row.names = F)
  write.csv(info,"sample.info.csv",row.names = F)
}
