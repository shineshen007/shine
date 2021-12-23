#' @title Column_Rows
#' @description print cols and rows
#' @param dir data path
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
Column_Rows <- function(dir='/Users/data.csv'){
  data.table::fread(dir) %>% dim(.)

}
