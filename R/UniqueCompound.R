#'unique compound
#'@export

UniqueCompound <- function(){
  data <- read.csv("data pathway.csv")
  n<-as.character(data$compound.name)
  name<-data.frame(unique(unlist(strsplit(n,";"))))
  colnames(name)<-"compound name"
  write.csv(name," unique metabolites.csv",row.names = F)
}
