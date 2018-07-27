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
  require(data.table)
  data <- fread("data_after_pre.csv")
  info <- fread("sample.info.csv")
  info<-setDF(info)
  data<- setDF(data)
  data<-data[,-1]
  info<-info[,-c(2:4)]
  write.csv(data,"data for dna.csv",row.names = F)
  write.csv(info,"info.csv",row.names = F)
}
