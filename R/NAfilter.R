#' NA filter
#' @export
#' @param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns.
NAfilter <- function(data){
  filter <- data[!is.na(data$compound.name),]
  write.csv(filter, "filter.NA.csv", row.names = FALSE)
}
