#' @title BindData
#' @description a function to combind data.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  bind pos and neg data.
#' @export
BindData <- function(){#bind pos and neg mode data

  data_pos <- data.table::fread("data pos.csv")
  data_pos<- data.table::setDF(data_pos)
  ID_pos <- data.table::fread("ID POS.csv")
  ID_pos <- data.table::setDF(ID_pos)
  data_neg <- data.table::fread("data neg.csv")
  data_neg<- data.table::setDF(data_neg)
  ID_neg <- data.table::fread("ID NEG.csv")
  ID_neg <- data.table::setDF(ID_neg)
  #r
  data_ID_pos<-cbind(ID_pos,data_pos)
  write.csv(data_ID_pos,"data_ID_pos.csv",row.names = F)
  data_ID_neg<-cbind(ID_neg,data_neg)
  write.csv(data_ID_neg,"data_ID_neg.csv",row.names = F)
  #
  ID_both<-rbind(ID_pos,ID_neg)
  write.csv(ID_both,"ID both.csv",row.names = F)
  #
  data_both<-rbind(data_ID_pos,data_ID_neg)
  write.csv(data_both,"data pos and neg.csv",row.names = F)
  data<-data_both[,-c(15:17)]
  write.csv(data,"data for next.csv",row.names = F)
}
