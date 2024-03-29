#' @title Metabolites_Classify
#' @description a function Metabolites_Classify.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param specias_pathway_database the specia of sample
#' @param specias_compound_database the specia of sample
#' @param data_name the data_name
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @import ggplot2
Metabolites_Classify <- function(data_name="data for classify.csv",
  specias_pathway_database= c(hsa.kegg.pathway,mmu.kegg.pathway),
  specias_compound_database = c(hsa_compound_ID,mmu_compound_ID)
){
  cat("Import data...\n")
  data <- data.table::fread(data_name)
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
  idi <- which(colnames(dalls)=="isotope")
  ida <- which(colnames(dalls)=="adduct")
  idg <- which(colnames(dalls)=="confidence")
  dcsore <- dalls[,-c(idn,idc,idi,ida,idg)]

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
  dsi <- add_column(dcsore,dall[,idns],.after = 'name')
  colnames(dsi)[2] <- "ID"
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
  dsic <- add_column(dsi,dallc[,idnc],.after = 'ID')
  colnames(dsic)[3] <- "compound.name"
  #
  ##unique isotope
  cat("unique the isotope ...\n")
  nci <- data$isotope
  cni <- strsplit(nci,split=';')
  a1 <- data.frame(NULL)
  a2 <- data.frame(NULL)
  a3 <- data.frame(NULL)
  for (i in 1:nr) {
    cat(i); cat(" ")
    if(length(cni[[i]]) != 1){
      rsci <- data[i,]
      times <- length(cni[[i]])
      for (j in 1:times) {
        rsci$isotope <- cni[[i]][j]
        a1 <- rbind(a1,rsci)
      }
    }else{
      a2 <- data[i,]
      a3 <- rbind(a3,a2)
    }
  }
  dalli <- rbind(a3,a1)
  idni <- which(colnames(dalli)=="isotope")
  dsici <- add_column(dsic,dalli[,idni],.after = 'compound.name')
  colnames(dsici)[4] <- "isotope"
  #
  ##unique adduct
  cat("unique the adduct ...\n")
  nca <- data$adduct
  cna <- strsplit(nca,split=';')
  a1 <- data.frame(NULL)
  a2 <- data.frame(NULL)
  a3 <- data.frame(NULL)
  for (i in 1:nr) {
    cat(i); cat(" ")
    if(length(cna[[i]]) != 1){
      rsca <- data[i,]
      times <- length(cna[[i]])
      for (j in 1:times) {
        rsca$adduct <- cna[[i]][j]
        a1 <- rbind(a1,rsca)
      }
    }else{
      a2 <- data[i,]
      a3 <- rbind(a3,a2)
    }
  }
  dalla <- rbind(a3,a1)
  idna <- which(colnames(dalla)=="adduct")
  dsica <- add_column(dsici,dalla[,idna],.after = 'isotope')
  colnames(dsica)[5] <- "adduct"
  ##unique confidence

  cat("unique the confidence ...\n")
  ncg <- data$confidence
  cng <- strsplit(ncg,split=';')
  a1 <- data.frame(NULL)
  a2 <- data.frame(NULL)
  a3 <- data.frame(NULL)
  for (i in 1:nr) {
    cat(i); cat(" ")
    if(length(cng[[i]]) != 1){
      rscg <- data[i,]
      times <- length(cng[[i]])
      for (j in 1:times) {
        rscg$confidence <- cng[[i]][j]
        a1 <- rbind(a1,rscg)
      }
    }else{
      a2 <- data[i,]
      a3 <- rbind(a3,a2)
    }
  }
  dallg <- rbind(a3,a1)
  idng <- which(colnames(dallg)=="confidence")
  dsicg <- add_column(dsica,dallg[,idng],.after = 'adduct')
  colnames(dsicg)[6] <- "confidence"
  #remove the isotope
  temp<- dsicg[c(grep("\\[M\\]",dsicg$isotope),
                 which(dsicg$isotope == "")),]
  #remove the score<.4
  idx <- which(temp$score < 0.4)
  df <- temp[-idx,]

  #unique compound
  qcinx <- grep('qc',colnames(df))
  qc <- df[,qcinx]
  qcm <- apply(qc, 1, mean)
  am <- add_column(df,qcm,.after = 'rt')
  max_uniq = aggregate(am[,"qcm"],list(am[,"compound.name"]),max,drop = FALSE)
  ic <- intersect(max_uniq$Group.1,am$compound.name)
  da <- am[match(ic,am$compound.name),]
  write.csv(da,'unique_compound.csv',row.names = F)
  #write.csv(df,'data after classify.csv',row.names = F)
  cat("classify the metabolites ...\n")
  #classify the metabolites
  mid <- as.character(da$ID)
  #mid <- as.character(mid)
  uid<-mid
  metabolite.id <- uid[which(uid %in% unique(unlist(specias_pathway_database)))]#filter the metabolites not in kegg

  #Allid <- as.character(as.matrix(metabolite.id))

  IDinPathway <- unlist(lapply(specias_pathway_database, function(x) {
    intersect(x, metabolite.id)#get the mapped metabolites number in each pathway
  }))
  # aID <- as.data.frame(IDinPathway[!duplicated(IDinPathway)])
  # i <- intersect(aID$`IDinPathway[!duplicated(IDinPathway)]`,metabolite.id)
  dt <- specias_compound_database[match(IDinPathway,specias_compound_database$id),]
  aa <- cbind(IDinPathway,dt)
  pmid <- aa[,-c(2:3)]
  colnames(pmid) <- c('ID','compound name')
  ######
  as <- da[,c('name','ID','compound.name','fc','p')]
  ass <- as[match(intersect(pmid$ID,as$ID),as$ID),]
  asp <- pmid[match(intersect(pmid$ID,ass$ID),pmid$ID),]
  aspd <- cbind(asp,ass)
  #######
  write.csv(aspd,'classfied metabolites.csv')
  ##Metabolite_Distribution_plot
  bk <- read.csv("classfied metabolites.csv")
  cda <- as.character(bk$X)
  cds <- unlist(strsplit(cda,";"))
  bboy <- grep('hsa',cds)#delete the rows hsa000101 etc
  bg <- cds[-bboy]
  bg <- as.data.frame(bg)
  ax <- plyr::ddply(bg, plyr::.(bg),summarize,number=length(bg))
  pnum <- nrow(unique(bg))
  ggplot2::ggplot(ax,ggplot2::aes(reorder(bg,number),number))+
    ggplot2::geom_bar(ggplot2::aes(fill=factor(1:pnum)),stat = "identity",position="dodge",width=0.8)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 1, vjust = .5))+
    ggplot2::guides(fill=FALSE)+
    ggplot2::coord_flip()+
    ggplot2::xlab('Pathway Name')+
    ggplot2::ylab(('The number of metabolites in each pathway'))+
    ggplot2::annotate("text", label = nrow(unique(bg)), x = 2, y = 20, size = 6)+
    ggplot2::annotate("text", label = 'pathway number', x = 5, y = 20, size = 6)+
    ggplot2::annotate("text", label = nrow(unique(bk)), x = 8, y = 20, size = 6)+
    ggplot2::annotate("text", label = 'metabolites number', x = 12, y = 20, size = 6)+
    ggplot2::geom_text(mapping = ggplot2::aes(label = ax$number),size=3,vjust=0.5,position = ggplot2::position_stack(vjust = 0.5))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),#remove ggplot2 background
                   panel.background = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"),legend.position = "none")+
    ggplot2::ggsave("Metabolite_Distribution_plot.png", width = 12, height = 8)
}
