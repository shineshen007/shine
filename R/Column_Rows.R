#' @title Column_Rows
#' @description print cols and rows
#' @param dir data path
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @export
Column_Rows <- function(dir='/Users/data.csv'){
  data.table::fread(dir) %>% dim(.)

}
