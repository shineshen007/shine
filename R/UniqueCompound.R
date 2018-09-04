#'unique compound
#'@export
#'@param data a dataframe must include compound name and ID.
UniqueCompound <- function(data){
  n<-as.character(data$compound.name)
  name<-unique(unlist(strsplit(n,";")))
  i<-as.character(data$ID)
  id<-unique(unlist(strsplit(i,";")))
  In<-cbind(id,name)
  colnames(In)<-c("kegg ID","compound name")
  write.csv(In," unique metabolites.csv",row.names = F)
}
