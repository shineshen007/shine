#' @title iMatch
#' @description a function to mtach data.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  get data you want.
#' @param name the column which you want to match
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
iMatch <- function(name = "name"){
   data1<-read.csv('data1.csv')
   data2<-read.csv('data2.csv')
   i <- intersect(data1$name,data2$name)
   dt <- data2[match(i,data2$name),]
   write.csv(dt,"results.csv",row.names = F)
 }















