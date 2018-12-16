#' @title IsotopeFilter
#' @description a function can filter isotope.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
IsotopeFilter <- function() {

  cat("Import data...\n")
  data <- data.table::fread("data.csv")
  data <- data.table::setDF(data)

  cat("Isotope filtering...\n")
  ###remove [M+n],\\ make [] lose the ability of function，
  isotope_filter<-function(data){
    temp<- data[c(grep("\\[M\\]",data$isotope),
                  which(data$isotope == "")),]
  }
  filter.isotope.data <-isotope_filter(data)
  write.csv(filter.isotope.data,"filter.isotope.csv",row.names = FALSE)

}
