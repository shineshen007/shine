#'unique compound
#'@export
#'@param data a dataframe must include compound name and ID.
UniqueCompound <- function(data){
  n<-as.character(data$compound.name)
  name<-unique(unlist(strsplit(n,";")))
  colnames(name)<-"compound name"
  write.csv(name," unique metabolites.csv",row.names = F)
}
