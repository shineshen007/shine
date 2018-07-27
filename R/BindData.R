#' @title BindData
#' @description a function to combind data.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
BindData <- function(){
  require(data.table)
  data_pos <- fread("data pos.csv")
  data_pos<- setDF(data_pos)
  MRN_pos <- fread("MRN POS.csv")
  MRN_pos <- setDF(MRN_pos)
  data_neg <- fread("data neg.csv")
  data_neg<- setDF(data_neg)
  MRN_neg <- fread("MRN NEG.csv")
  MRN_neg <- setDF(MRN_neg)
  #
  data_MRN_pos<-cbind(MRN_pos,data_pos)
  write.csv(data_MRN_pos,"data_MRN_pos.csv",row.names = F)
  data_MRN_neg<-cbind(MRN_neg,data_neg)
  write.csv(data_MRN_neg,"data_MRN_neg.csv",row.names = F)
  #
  data_both<-rbind(data_MRN_pos,data_MRN_neg)
  write.csv(data_both,"data pos and neg.csv",row.names = F)
  data<-data_both[,-c(4:6,11:17)]
  write.csv(data,"data for next.csv",row.names = F)
}
