#' @title DataFilter
#' @description a function to filter p value more than 0.05 or any value you can choose,
#' foldchange less than 1.2 and more than 0.8 or any value you can choose,
#' and vip less than 1 or any value you can choose
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param p p value
#' @param comp1 component 1 in PLSDA
#' @param comp2 component 2 in PLSDA
#' @param fc_toplimit top limit of fold change
#' @param fc_lowerlimit lower limit of fold change
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
DataFilter <- function(p=0.05,comp1=1,comp2=1,
                        fc_toplimit=4/3,fc_lowerlimit=3/4){
  require(data.table)
  data <- fread("data_pfc_vip.csv")
  data <- setDF(data)
        idx.vip1<-which(data$`comp 1`<comp1)
        data<-data[-idx.vip1,]
        idx.vip2<-which(data$`comp 2`<comp2)
        data<-data[-idx.vip2,]
        idx.fc1<-which(data$fc<fc_lowerlimit)
        data1<-data[idx.fc1,]
        idx.fc2<-which(data$fc>fc_toplimit)
        data2<-data[idx.fc2,]
        data_mid<-rbind(data1,data2)
        idx.p<-which(data_mid$p>p)
        data_final<-data_mid[-idx.p,]
        write.csv(data_final,"data final.csv",row.names = FALSE)
  }