#' @title PathwayEnrich
#' @description a function can match differetiate metabolites in differentiate pathway.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param specias_pathway_database the specia of sample
#' @param specias_compound_database the specia of sample
#' @param keggmap generate file for keggmap
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
PathwayEnrich <- function(specias_pathway_database= c(hsa.kegg.pathway,mmu.kegg.pathway),
                          specias_compound_database = c(hsa_compound_ID,mmu_compound_ID),
                          keggmap = TRUE
){
  data <- read.csv('data pathway.csv')

  mid <- data$ID
  mid <- as.character(mid)
  uid<-unique(unlist(strsplit(mid,";")))
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
    nsd <- match(sd,aID$`IDinPathway[!duplicated(IDinPathway)]`)
    aa <- cbind(aID[-nsd,],dt)
  }else(aa <- cbind(aID,dt))

if(keggmap){
  ###
  nr <- nrow(data)
  nid <- data$ID
  cnid <- strsplit(nid,split=';')
  a1 <- data.frame(NULL)
  a2 <- data.frame(NULL)
  a3 <- data.frame(NULL)
  for (i in 1:nr) {
    cat(i); cat(" ")
    if(length(cnid[[i]]) != 1){
      rss <- data[i,]
      times <- length(cnid[[i]])
      for (j in 1:times) {
        rss$ID <- cnid[[i]][j]
        a1 <- rbind(a1,rss)
      }
    }else{
      a2 <- data[i,]
      a3 <- rbind(a3,a2)
    }
  }
  dall <- rbind(a3,a1)
  idx <- intersect(aa$id,dall$ID)
  alld <- dall[match(idx,dall$ID),]
  aaldd <- cbind(aa,alld)
  rn<-nrow(aaldd)
  lfc<-which(aaldd$fc<1)
  tfc<-which(aaldd$fc>1)
  ab<-data.frame(1:rn)
  colnames(ab)<-"colour"
  ab[lfc,]<-"green"
  ab[tfc,]<-"red"
  dfc<-cbind(ab,aaldd)
  dfcc <- dfc[,c(2,1)]
  ##create a folder for analysis
  dir.create('PathwayEnrich')
  setwd('PathwayEnrich')
  write.csv(dfc,'differentiate metabolites in pathway.csv')
  write.csv(dfcc,'metabolites map to pathway.csv')

}

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
  w <- which(dat$p.value>0.05)
  row = w[1]+5
  dat <- dat[1:row,]
  group <- ifelse(dat$p < 0.05,"sig", "not sig")
  pb<- ggplot2::ggplot(dat,ggplot2::aes(reorder(pathway,-p),-log10(p)))+##-p control the order
    ggplot2::geom_bar(aes(fill=group),stat = "identity",position="dodge",width=0.8)+
    ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
                   panel.background = element_blank(),axis.line = element_line(colour = "black"),
                   legend.position = "none",axis.text.y = element_text(size = 14),#the font size of axis
                   axis.title.x = element_text(size = 14),#the font size of axis title
                   axis.title.y = element_text(size = 14))+
    ggplot2::geom_hline(yintercept= 1.30103,#draw a line at p=0.05
                        lty=4,col="grey21",lwd=1)+
    ggplot2::labs(title="Pathway Enrichment Analysis")+
    ggplot2::coord_flip()+
    ggplot2::xlab('Pathway')+
    ggplot2::ggsave("PathwayBarplot.png", width = 12, height = 8)
  export::graph2ppt(x=pb,file='PathwayBarplot.pptx',width=12,height=8)
  ##back origin work directory
  setwd("..//")
}

