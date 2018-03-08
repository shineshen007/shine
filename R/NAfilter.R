#### NA filter
NA.filter <- function(data){
  filter <- data[!is.na(data$compound.name),]
  write.csv(filter, "filter.NA.csv", row.names = FALSE)
}
