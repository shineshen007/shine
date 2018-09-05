#' @title PathwayMatch
#' @description a function can match differetiate metabolites in differentiate pathway.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param group group set.
#' @param hsa hsa or mmu
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' library(Shine)
#' PathwayMatch(group = c("S","P"))
#' }
PathwayMatch<-function(group=c("case","control"),hsa=TRUE){
  require(xlsx)
  data <- read.csv("Quantitative.pathway.metabolite.result.csv")
  sample.info <- read.csv("sample.info.csv")
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]

  sample<-data[,match(sample.name,colnames(data))]

  name <- as.character(data[,"compound.name"])

  class<- sample.info[,"group"]
  pathway.p<-read.csv("Pathway.enrichment.analysis.csv")
  nume<-length(which(pathway.p$FDR<0.05))
  path1<-as.character(pathway.p[1:nume,1])
  c<-data.frame(NULL)
  pathwayy<-function(){
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

    p <- p.adjust(p = p, method = "fdr",n=length(p))

    f<-as.data.frame(fc)
    ##get kegg code
    dfc<-data[,1:2]
    dfc<-cbind(dfc,f)
    rn<-dim(dfc)[1]
    lfc<-which(dfc$fc<1)
    tfc<-which(dfc$fc>1)
    ab<-data.frame(1:rn)
    colnames(ab)<-"colour"
    ab[lfc,]<-"green"
    ab[tfc,]<-"red"
    dfc<-cbind(dfc,ab)
    dfc<-dfc[,c(1,3,2,4)]
    write.csv(dfc,"kegg.csv",row.names = F)
    #
    pvalue<-as.data.frame(p)
    data_vol<-cbind(name,f,pvalue)
    #get the vol.csc file ,merge with Quantitative.pathway.metabolite.result.csv
    data.path<-cbind(data_vol,data)
    y<-data.path$ID
    int<-intersect(y,pathway)
    n<-as.data.frame(int)
    num<-nrow(n)
    match<-data.frame(data.path$name[match(int,y)])
    rowname<-as.data.frame(c(rep(c(pathway.name),num)))
    path<-cbind(rowname,match)
    match<-data.frame(data.path$fc[match(int,y)])
    pathway.pvalue<-match(pathway.name,pathway.p$pathway)
    pvalue<-pathway.p[pathway.pvalue,4]
    pp<-as.data.frame(c(rep(c(pvalue),num)))
    cyto<-cbind(pp,path,match)
    colnames(cyto)<-c("pathway.p","pathway","metabolites","fc")
    cyto<-cyto
  }
  if(hsa){
  for (i in 1:nume){
    pathway.name<-path1[i]
    pathway.d<-hsa.kegg.pathway[pathway.name]
    pathway<-pathway.d[[1]]
    a<-NULL
    d<-pathwayy()
    temp.c <- rbind(a,d)
    c <- rbind(c, temp.c)
  }
  }else{
    for (i in 1:nume){
      pathway.name<-path1[i]
      pathway.d<-mmu.kegg.pathway[pathway.name]
      pathway<-pathway.d[[1]]
      a<-NULL
      d<-pathwayy()
      temp.c <- rbind(a,d)
      c <- rbind(c, temp.c)
}
  }
  write.xlsx(c,"cytoscape.xlsx",row.names = F)

}
