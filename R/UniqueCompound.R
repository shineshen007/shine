#'unique compound
#'@export

UniqueCompound <- function(data_name="data pathway.csv"){
  data <- read.csv(data_name)
  n<-as.character(data$compound.name)
  name<-data.frame(unique(unlist(strsplit(n,";"))))
  colnames(name)<-"compound name"
  write.csv(name,"unique metabolites.csv",row.names = F)
}
