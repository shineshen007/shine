#' @title IsotopeFilter
#' @description a function can filter isotope.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
IsotopeFilter <- function() {

  require(data.table)
  cat("Import data...\n")
  data <- fread("data.csv")
  data <- setDF(data)

  cat("Isotope filtering...\n")
  ###remove [M+n],\\ make [] lose the ability of function，
  isotope_filter<-function(data){
    temp<- data[c(grep("\\[M\\]",data$isotopes),
                  which(data$isotopes == "")),]
  }
  filter.isotope.data <-isotope_filter(data)
  write.csv(filter.isotope.data,"filter.isotope.csv",row.names = FALSE)

}
