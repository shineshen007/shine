### unique compound
unique_compound<-function(data){
  i<-as.character(data$compound.name)
  a<-unique(unlist(strsplit(i,";")))
  write.csv(a," unique metabolites.csv",row.names = F)
}
