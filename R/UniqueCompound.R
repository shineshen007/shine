#'unique compound
#'@export
#'@param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns.
UniqueCompound<-function(data){
  i<-as.character(data$compound.name)
  a<-unique(unlist(strsplit(i,";")))
  write.csv(a," unique metabolites.csv",row.names = F)
}
