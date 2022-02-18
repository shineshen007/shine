#' @title iMatch
#' @description a function to mtach data.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  get data you want.
#' @param d1 the position of the file in dir()
#' @param d2 the position of the file in dir()
#' @param name the column which you want to match
#' @export
iMatch <- function(d1=1,
                   d2=2,#the data you want to extract
                   name = "name"#the column you want to match
                   ){
   data1<-read.csv(dir()[d1])
   data2<-read.csv(dir()[d2])
   id <- which(colnames(data1)==name)
   i <- intersect(data1[,id],data2[,id])
   dt <- data2[match(i,data2[,id]),]
   write.csv(dt,"results.csv",row.names = F)
}















