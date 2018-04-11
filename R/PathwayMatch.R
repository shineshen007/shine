#' @title PathwayMatch
#' @description a function can match differetiate metabolites in differentiate pathway.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data a dataframe from MetDNA named Quantitative.pathway.metabolite.result.
#' @param sample.info a dataframe include sample.name,injection.order,
#' class,batch and group columns.
#' @param group group set.
#' @param pcorrect default is TRUE.
#' @param pathway.name it is the pathway you want to match.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' library(Shine)
#' data("hsa.kegg.pathway",package="Shine")
#' data("mmu.kegg.pathway",package="Shine")
#' load(h)
#' load(m)
#' pathway<-hsa.kegg.pathway$'Cysteine and methionine metabolism'
#' PathwayMatch(data = data,sample.info = sample.info,group = c("S","P"),
#' pcorrect = F,pathway.name = "Cysteine and methionine metabolism")
#' }
PathwayMatch <- function(data = NULL,sample.info = NULL,
                       group = c("case","control"),
                       pcorrect = TRUE,pathway.name = pathway){
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]

  sample<-data[,match(sample.name,colnames(data))]

  name <- as.character(data[,"compound.name"])

  class<- sample.info[,"group"]

  group1.index <- which(class == group[1])
  group2.index <- which(class == group[2])

  fc <- apply(sample,1,function(x) {
    median(x[group2.index]+0.1)/ median(x[group1.index]+0.1)
  })

  t.test <- apply(sample, 1, function(x) {
    t.test(x[group1.index], x[group2.index])
  })

  p <- unlist(lapply(t.test, function(x)
    x$p.value))
  if(pcorrect){
    p <- p.adjust(p = p, method = "fdr",n=length(p))
  }
  f<-as.data.frame(fc)
  pvalue<-as.data.frame(p)
  data_vol<-cbind(name,f,pvalue)
  #get the vol.csc file ,merge with Quantitative.pathway.metabolite.result.csv
  data.path<-cbind(data_vol,data)
  y<-data.path$ID
  i<-intersect(y,pathway)
  n<-as.data.frame(i)
  num<-nrow(n)
  match<-data.frame(data.path$name[match(i,y)])
  rowname<-as.data.frame(c(rep(c(pathway.name),num)))
  path<-cbind(rowname,match)
  match<-data.frame(data.path$fc[match(i,y)])
  cyto<-cbind(path,match)
  colnames(cyto)<-c("pathway","metabolites","fc")
  write.csv(cyto,"cytoscape.csv",row.names = F)
}
