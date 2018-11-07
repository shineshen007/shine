#' @title BindData
#' @description a function to combind data.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  bind pos and neg data.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
BindData <- function(){#bind pos and neg mode data
  require(data.table)
  data_pos <- fread("data pos.csv")
  data_pos<- setDF(data_pos)
  ID_pos <- fread("ID POS.csv")
  ID_pos <- setDF(ID_pos)
  data_neg <- fread("data neg.csv")
  data_neg<- setDF(data_neg)
  ID_neg <- fread("ID NEG.csv")
  ID_neg <- setDF(ID_neg)
  #
  data_ID_pos<-cbind(ID_pos,data_pos)
  write.csv(data_ID_pos,"data_ID_pos.csv",row.names = F)
  data_ID_neg<-cbind(ID_neg,data_neg)
  write.csv(data_ID_neg,"data_ID_neg.csv",row.names = F)
  #
  data_both<-rbind(data_ID_pos,data_ID_neg)
  write.csv(data_both,"data pos and neg.csv",row.names = F)
  data<-data_both[,-c(4:6,11:17)]
  write.csv(data,"data for next.csv",row.names = F)
}
