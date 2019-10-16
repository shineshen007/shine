#' @title PathwayEnrich
#' @description a function can match differetiate metabolites in differentiate pathway.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param specias_pathway_database the specia of sample
#' @param specias_compound_database the specia of sample
#' @param keggmap default is TRUE
#' @param group group info.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
PathwayEnrich <- function(specias_pathway_database= c(hsa.kegg.pathway,mmu.kegg.pathway),
                          specias_compound_database = c(hsa_compound_ID,mmu_compound_ID),
                          keggmap = TRUE,
                          group=c("case","control")
){
  da <- fread('data after classify.csv')
  da<-setDF(da)
  #unique ID
  max_score = aggregate(da[,"score"],list(da[,"compound.name"]),max,drop = FALSE)
  ic <- intersect(max_score$x,da$score)
  data <- da[match(ic,da$score),]
  sample.info <- read.csv("sample.info.csv")
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  qc.name<-sample.info$sample.name[sample.info$class=="QC"]

  sample<-data[,match(sample.name,colnames(data))]
  qc<-data[,match(qc.name,colnames(data))]

  data_pfc<- as.matrix(cbind(qc,sample))


  name <- as.character(data[,"compound.name"])

  class<- sample.info[,"group"]

  group1.index <- which(class == group[1])
  group2.index <- which(class == group[2])

  fc <- apply(data_pfc,1,function(x) {
    median(x[group1.index]+0.1)/ median(x[group2.index]+0.1)
  })

  ptest <- apply(data_pfc, 1, function(x) {
    wilcox.test(x[group1.index], x[group2.index])
  })

  p <- unlist(lapply(ptest, function(x)
    x$p.value))

  p <- p.adjust(p = p, method = "fdr",n=length(p))

  f<-as.data.frame(fc)

  #
  pvalue<-as.data.frame(p)
  data_vol<-cbind(name,f,pvalue)
  #get the vol.csc file ,merge with Quantitative.pathway.metabolite.result.csv
  data.path<-cbind(data_vol,data)
  mp <- which(data.path$p>0.05)
  data.path <- data.path[-mp,]
  ##create a folder for analysis
  dir.create('PathwayEnrich')
  setwd('PathwayEnrich')

  if(keggmap){
    rn<-nrow(data.path)
    lfc<-which(data.path$fc<1)
    tfc<-which(data.path$fc>1)
    ab<-data.frame(1:rn)
    colnames(ab)<-"colour"
    ab[lfc,]<-"green"
    ab[tfc,]<-"red"
    dfc<-cbind(data.path[,"ID"],ab)

    write.csv(dfc,'metabolites map to pathway.csv')
  }
  #pathwayenrich
  mid <- data.path$ID
  # mid <- as.character(mid)
  uid<-mid
  metabolite.id <- uid[which(uid %in% unique(unlist(specias_pathway_database)))]#filter the metabolites not in specia

  ALLm <- unname(unique(unlist(specias_pathway_database)))
  ALLm <- as.character(as.matrix(ALLm))
  SIGm <- as.character(as.matrix(metabolite.id))
  num_all <- length(ALLm)
  num_sig <- length(SIGm)

  pall <- unlist(lapply(specias_pathway_database, function(x) {
    length(intersect(x, ALLm))#get the metabolites number in each pathway
  }))

  psig <- unlist(lapply(specias_pathway_database, function(x) {
    length(intersect(x, SIGm))#get the mapped metabolites number in each pathway
  }))

  IDinPathway <- unlist(lapply(specias_pathway_database, function(x) {
    intersect(x, SIGm)#get the mapped metabolites number in each pathway
  }))
  aID <- as.data.frame(IDinPathway[!duplicated(IDinPathway)])
  i <- intersect(aID$`IDinPathway[!duplicated(IDinPathway)]`,specias_compound_database$id)
  dt <- specias_compound_database[match(i,specias_compound_database$id),]
  sd <- setdiff(aID$`IDinPathway[!duplicated(IDinPathway)]`,dt$id)
  if(length(sd) != 0){
    nsd <- which(aID$`IDinPathway[!duplicated(IDinPathway)]`%in%sd)
    aa <- cbind(aID[-nsd,],dt)
  }else(aa <- cbind(aID,dt))

  ###

  P<-NaN
  for (i in 1:length(pall)){
    ball <- psig[i]
    P[i] <- phyper(q = ball, m = pall[i], n = num_all - pall[i],
                   k = num_sig, lower.tail = FALSE)
  }
  BHP <- p.adjust(P, method="BH")
  FDR <- p.adjust(P, method="fdr")

  PBF <- cbind(P,BHP, FDR)
  rownames(PBF) <- colnames(psig)
  colnames(PBF) <- c("p.value","q.value", "FDR")

  ##calculate the impact of pathway
  info <- lapply(specias_pathway_database, function(module) {
    overlap.number <- length(intersect(module, metabolite.id))
    pathway.number <- length(module)
    c(pathway.number, overlap.number)
  })

  info <- do.call(rbind, info)
  colnames(info) <- c("Pathway.length", "Overlap")

  info <- data.frame(PBF, info, stringsAsFactors = FALSE)

  info <- info[order(info[,1]),]

  write.csv(info,'pathway impact.csv')
  #draw barplot
  dat <- read.csv('pathway impact.csv')

  colnames(dat) <- c('pathway', "p.value","q.value","p","Pathway.length","Overlap" )
  w <- which(dat$p>0.05)
  row = w[1]+5
  dat <- dat[1:row,]
  pa <- as.character(dat$pathway)
  pid<-unique(unlist(strsplit(pa,";")))
  an <- seq(2,length(pid),2)
  pid <- data.frame(pid[-an])
  colnames(pid) <- 'path'
  datt <- cbind(pid,dat)
  group <- ifelse(dat$p < 0.05,"sig", "not sig")
  cat("Draw pathway barplot...\n")
  pb<- ggplot2::ggplot(datt,ggplot2::aes(reorder(path,-p),-log10(p)))+##-p control the order
    ggplot2::geom_bar(ggplot2::aes(fill=group),stat = "identity",position="dodge",width=0.8)+
    ggplot2::theme(panel.grid.major =ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),#remove ggplot2 background
                   panel.background = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"),
                   legend.position = "none",axis.text.y = ggplot2::element_text(size = 16),#the font size of axis
                   axis.title.x = ggplot2::element_text(size = 18),#the font size of axis title
                   axis.title.y = ggplot2::element_text(size = 18))+
    ggplot2::geom_hline(yintercept= 1.30103,#draw a line at p=0.05
                        lty=4,col="grey21",lwd=1)+
    ggplot2::labs(title="Pathway Enrichment Analysis")+
    ggplot2::coord_flip()+
    ggplot2::xlab('Pathway')+
    ggplot2::ggsave("PathwayBarplot.png", width = 12, height = 8)
  export::graph2ppt(x=pb,file='PathwayBarplot.pptx',height=7,width=9)
  cat("Generate cytoscape data...\n")
  #cytoscape
  pathway.p<-dat
  nume<-length(pathway.p$p<0.05)
  path1<-as.character(pathway.p[1:nume,1])
  a<-NULL
  ct<-data.frame(NULL)
  for (i in 1:nume){
    #get the pathway
    pathway.name<-path1[i]
    pathway.d<-specias_pathway_database[pathway.name]
    #get the ID in pathway
    pathway<-pathway.d[[1]]
    ix <- intersect(pathway,data.path$ID)
    pathway_name<-as.data.frame(c(rep(c(pathway.name),length(ix))))
    pathway_p<-as.data.frame(c(rep(c(pathway.p[i,'p']),length(ix))))
    foldchange <- data.path$fc[match(ix,data.path$ID)]
    metabolite <- data.path$compound.name[match(ix,data.path$ID)]
    cyto <- cbind(pathway_p,pathway_name,metabolite,foldchange)

    ct <- rbind(ct, cyto)

  }
  colnames(ct)<-c("pathway.p","pathway","metabolites","foldchange")
  xlsx::write.xlsx(ct,"cytoscape.xlsx",row.names = F)

  ##back origin work directory
  setwd("..//")
}

