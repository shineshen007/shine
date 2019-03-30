#' @title Metabolites_Classify
#' @description a function Metabolites_Classify.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param specias_pathway_database the specia of sample
#' @param specias_compound_database the specia of sample
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @import ggplot2
#' @examples
#' \donttest{
#' }
Metabolites_Classify <- function(
  specias_pathway_database= c(hsa.kegg.pathway,mmu.kegg.pathway),
  specias_compound_database = c(hsa_compound_ID,mmu_compound_ID)
){
  cat("Import data...\n")
  data <- data.table::fread("data for classify.csv")
  data <- data.table::setDF(data)
  cat("filter the metabolites' score lower than 0.4...\n")
  ##filter the score lower than 0.4
  nr <- nrow(data)
  nscore <- data$score
  cs <- strsplit(nscore,split=';')
  a1 <- data.frame(NULL)
  a2 <- data.frame(NULL)
  a3 <- data.frame(NULL)
  for (i in 1:nr) {
    cat(i); cat(" ")
    if(length(cs[[i]]) != 1){
      rss <- data[i,]
      times <- length(cs[[i]])
      for (j in 1:times) {
        rss$score <- cs[[i]][j]
        a1 <- rbind(a1,rss)
      }
    }else{
      a2 <- data[i,]
      a3 <- rbind(a3,a2)
    }
  }
  dalls <- rbind(a3,a1)
  idn <- which(colnames(dalls)=="ID")
  idc <- which(colnames(dalls)=="compound.name")
  dcsore <- dalls[,-c(idn,idc)]
  ##unique ID
  cat("unique the ID ...\n")
  nid <- data$ID
  cnid <- strsplit(nid,split=';')
  a1 <- data.frame(NULL)
  a2 <- data.frame(NULL)
  a3 <- data.frame(NULL)
  for (i in 1:nr) {
    cat(i); cat(" ")
    if(length(cnid[[i]]) != 1){
      rsi <- data[i,]
      times <- length(cnid[[i]])
      for (j in 1:times) {
        rsi$ID <- cnid[[i]][j]
        a1 <- rbind(a1,rsi)
      }
    }else{
      a2 <- data[i,]
      a3 <- rbind(a3,a2)
    }
  }
  dall <- rbind(a3,a1)
  idns <- which(colnames(dall)=="ID")
  dsi <- cbind(dall[,idns],dcsore)
  colnames(dsi)[1] <- "ID"
  ##unique compound
  cat("unique the compound ...\n")
  nc <- data$compound.name
  cnc <- strsplit(nc,split=';')
  a1 <- data.frame(NULL)
  a2 <- data.frame(NULL)
  a3 <- data.frame(NULL)
  for (i in 1:nr) {
    cat(i); cat(" ")
    if(length(cnc[[i]]) != 1){
      rsc <- data[i,]
      times <- length(cnc[[i]])
      for (j in 1:times) {
        rsc$compound.name <- cnc[[i]][j]
        a1 <- rbind(a1,rsc)
      }
    }else{
      a2 <- data[i,]
      a3 <- rbind(a3,a2)
    }
  }
  dallc <- rbind(a3,a1)
  idnc <- which(colnames(dallc)=="compound.name")
  dsic <- cbind(dallc[,idnc],dsi)
  colnames(dsic)[1] <- "compound.name"
  idx <- which(dsic$score < 0.4)
  data <- dsic[-idx,]
  write.csv(data,'data for classify.csv')
  cat("classify the metabolites ...\n")
  #classify the metabolites
  mid <- data$ID
  mid <- as.character(mid)
  uid<-unique(unlist(strsplit(mid,";")))
  metabolite.id <- uid[which(uid %in% unique(unlist(specias_pathway_database)))]#filter the metabolites not in kegg

  Allid <- as.character(as.matrix(metabolite.id))

  IDinPathway <- unlist(lapply(specias_pathway_database, function(x) {
    intersect(x, Allid)#get the mapped metabolites number in each pathway
  }))
  aID <- as.data.frame(IDinPathway[!duplicated(IDinPathway)])
  i <- intersect(aID$`IDinPathway[!duplicated(IDinPathway)]`,metabolite.id)
  dt <- specias_compound_database[match(i,specias_compound_database$id),]
  aa <- cbind(aID,dt)
  pmid <- aa[,-c(2:3)]
  colnames(pmid) <- c('ID','compound name')
  write.csv(pmid,'classfied metabolites.csv')
  ##Metabolite_Distribution_plot
  bk <- read.csv("classfied metabolites.csv")
  cda <- as.character(bk$X)
  cds <- unlist(strsplit(cda,";"))
  bboy <- grep('hsa',cds)
  bg <- cds[-bboy]
  bg <- as.data.frame(bg)
  ax <- plyr::ddply(bg, .(bg),summarize,number=length(bg))
  pnum <- nrow(unique(bg))
  ggplot2::ggplot(ax,aes(reorder(bg,number),number))+
    ggplot2::geom_bar(aes(fill=factor(1:pnum)),stat = "identity",position="dodge",width=0.8)+
    ggplot2::theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = .5))+
    ggplot2::guides(fill=FALSE)+
    ggplot2::coord_flip()+
    ggplot2::xlab('Pathway Name')+
    ggplot2::ylab(('The metabolites number in each pathway'))+
    ggplot2::annotate("text", label = nrow(unique(bg)), x = 2, y = 20, size = 6)+
    ggplot2::annotate("text", label = 'pathway number', x = 5, y = 20, size = 6)+
    ggplot2::annotate("text", label = nrow(unique(bk)), x = 8, y = 20, size = 6)+
    ggplot2::annotate("text", label = 'metabolites number', x = 12, y = 20, size = 6)+
    ggplot2::geom_text(mapping = aes(label = ax$number),size=3,vjust=0.5,position = position_stack(vjust = 0.5))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),#remove ggplot2 background
                   panel.background = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"),legend.position = "none")+
    ggplot2::ggsave("Metabolite_Distribution_plot.png", width = 12, height = 8)
}
