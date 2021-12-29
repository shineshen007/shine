#' @title PathwayEnrich
#' @description a function to do PathwayAnalysis
#' @author Xia Shen
#' \email{qq951633542@@163.com}
#' @param font_size 20
#' @param seed default is TRUE
#' @param All_ID default is FALSE
#' @param grade1_2 default is FALSE
#' @param specias_pathway_database the specia of sample
#' @param specias_compound_database the specia of sample
#' @param color color
#' @return  All the results can be got form other functions and instruction.
#' @export
PathwayEnrich <- function(specias_pathway_database= c(hsa.kegg.pathway,mmu.kegg.pathway),
                          specias_compound_database = c(hsa_compound_ID,mmu_compound_ID),
                          font_size = 20,
                          color = c("#3399FF","#FF3333"),
                          All_ID = FALSE,
                          seed = TRUE,
                          grade1_2 = FALSE
){
  pacman::p_load(ggplotify,readr,data.table,export,magrittr)
  if(All_ID){
    data <- read_csv('data_pathway.csv')
    if(!file.exists('PathwayEnrich_All_ID')){
      dir.create('PathwayEnrich_All_ID')
    }
    mid <- data$ID
    metabolite.id <- mid[which(mid %in% unique(unlist(specias_pathway_database)))]#filter the metabolites not in specia

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
    ###
    P<-NaN
    for (i in 1:length(pall)){
      ball <- psig[i]
      P[i] <- phyper(q = ball, m = pall[i], n = num_all - pall[i],
                     k = num_sig, lower.tail = FALSE)
    }
    BHP <- p.adjust(P, method="bonferroni")
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
    setwd('PathwayEnrich_All_ID')
    info <- info[order(info[,1]),]
    write.csv(info,'pathway_results.csv')
    #o('pathway_impact_76_pathway_hsa.all.ID.csv')
    sig_path <- dplyr::filter(info,q.value<0.05)
    ip <- lapply(specias_pathway_database, function(x){
      intersect(x,mid)
    })

    for (i in setdiff(names(specias_pathway_database),rownames(sig_path))) {
      ip[[i]]=NULL
    }

    for (i in 1:length(ip)) {
      ip[[i]] <- data[match(ip[[i]],data$ID),]
    }
    openxlsx::write.xlsx(ip,'Id_pathway.xlsx')
    #
    data <- read.csv("pathway_results.csv")
    colnames(data)[3] <- 'p'
    row=sum(data$p<0.05)+3
    data <- data[1:row,]
    group <- ifelse(data$p < 0.05,"sig", "not sig")
    pb <- ggplot2::ggplot(data,ggplot2::aes(reorder(X,-p),-log10(p)))+##-p control the order
      ggplot2::geom_bar(aes(fill=group),stat = "identity",position="dodge",width=0.8)+
      scale_fill_manual(values = color)+
      ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
                     panel.background = element_blank(),axis.line = element_line(colour = "black"),
                     legend.position = "none",axis.text.y = element_text(size = font_size),
                     axis.text.x = element_text(size = font_size),#the font size of axis
                     axis.title.x = element_text(size = font_size),#the font size of axis title
                     axis.title.y = element_text(size = font_size))+
      ggplot2::geom_hline(yintercept= 1.30103,#draw a line at p=0.05
                          lty=4,col="grey21",lwd=1)+
      ggplot2::labs(title="Pathway Enrichment Analysis")+
      ggplot2::coord_flip()+
      ggplot2::xlab('Pathway')
    save(pb,file = 'barplot.Rda')
    export::graph2ppt(x=pb,file='pathway.pptx',height=7,width=9)
    dev.off()
    setwd('..//')
  }
  #####
  if(seed){
    dir.create('PathwayEnrich_seed')
    #pathwayenrich
    data <- read_csv('data_pathway.csv')
    data <- dplyr::filter(data,Annotation.type=='seed')
    #seed <- dplyr::filter(da,confidence=='grade1'|confidence=='grade2')
    mid <- data$ID
    metabolite.id <- mid[which(mid %in% unique(unlist(specias_pathway_database)))]#filter the metabolites not in specia

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
    ###
    P<-NaN
    for (i in 1:length(pall)){
      ball <- psig[i]
      P[i] <- phyper(q = ball, m = pall[i], n = num_all - pall[i],
                     k = num_sig, lower.tail = FALSE)
    }
    BHP <- p.adjust(P, method="bonferroni")
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
    setwd('PathwayEnrich_seed')
    info <- info[order(info[,1]),]
    write.csv(info,'pathway_results.csv')
    #o('pathway_impact_76_pathway_hsa.all.ID.csv')
    ##id match to pathway
    #compound_db <- read_csv('/Users/shineshen/Downloads/demo/kegg/old/kegg_compound_unique.csv')
    sig_path <- dplyr::filter(info,q.value<0.05)
    ip <- lapply(specias_pathway_database, function(x){
      intersect(x,mid)
    })

    for (i in setdiff(names(specias_pathway_database),rownames(sig_path))) {
      ip[[i]]=NULL
    }

    for (i in 1:length(ip)) {
      ip[[i]] <- data[match(ip[[i]],data$ID),]
    }
    openxlsx::write.xlsx(ip,'Id_pathway.xlsx')
    #
    data <- read.csv("pathway_results.csv")
    colnames(data)[3] <- 'p'
    row=sum(data$p<0.05)+3
    data <- data[1:row,]
    group <- ifelse(data$p < 0.05,"sig", "not sig")
    pb <- ggplot2::ggplot(data,ggplot2::aes(reorder(X,-p),-log10(p)))+##-p control the order
      ggplot2::geom_bar(aes(fill=group),stat = "identity",position="dodge",width=0.8)+
      scale_fill_manual(values = color)+
      ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
                     panel.background = element_blank(),axis.line = element_line(colour = "black"),
                     legend.position = "none",axis.text.y = element_text(size = font_size),
                     axis.text.x = element_text(size = font_size),#the font size of axis
                     axis.title.x = element_text(size = font_size),#the font size of axis title
                     axis.title.y = element_text(size = font_size))+
      ggplot2::geom_hline(yintercept= 1.30103,#draw a line at p=0.05
                          lty=4,col="grey21",lwd=1)+
      ggplot2::labs(title="Pathway Enrichment Analysis")+
      ggplot2::coord_flip()+
      ggplot2::xlab('Pathway')
    save(pb,file = 'barplot.Rda')
    export::graph2ppt(x=pb,file='pathway.pptx',height=7,width=9)
    dev.off()
    setwd('..//')
  }
  ###########
  if(grade1_2){
    dir.create('PathwayEnrich_grade1_2')
    #pathwayenrich
    #data <- dplyr::filter(data,Annotation.type=='seed')
    data <- read_csv('data_pathway.csv')
    data <- dplyr::filter(data,confidence=='grade1'|confidence=='grade2')
    mid <- data$ID
    metabolite.id <- mid[which(mid %in% unique(unlist(specias_pathway_database)))]#filter the metabolites not in specia

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

    P<-NaN
    for (i in 1:length(pall)){
      ball <- psig[i]
      P[i] <- phyper(q = ball, m = pall[i], n = num_all - pall[i],
                     k = num_sig, lower.tail = FALSE)
    }
    BHP <- p.adjust(P, method="bonferroni")
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
    setwd('PathwayEnrich_grade1_2')
    info <- info[order(info[,1]),]
    write.csv(info,'pathway_results.csv')
    #o('pathway_impact_76_pathway_hsa.all.ID.csv')
    sig_path <- dplyr::filter(info,q.value<0.05)
    ip <- lapply(specias_pathway_database, function(x){
      intersect(x,mid)
    })

    for (i in setdiff(names(specias_pathway_database),rownames(sig_path))) {
      ip[[i]]=NULL
    }

    for (i in 1:length(ip)) {
      ip[[i]] <- data[match(ip[[i]],data$ID),]
    }
    openxlsx::write.xlsx(ip,'Id_pathway.xlsx')
    #
    data <- read.csv("pathway_results.csv")
    colnames(data)[3] <- 'p'
    row=sum(data$p<0.05)+3
    data <- data[1:row,]
    group <- ifelse(data$p < 0.05,"sig", "not sig")
    pb <- ggplot2::ggplot(data,ggplot2::aes(reorder(X,-p),-log10(p)))+##-p control the order
      ggplot2::geom_bar(aes(fill=group),stat = "identity",position="dodge",width=0.8)+
      scale_fill_manual(values = color)+
      ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
                     panel.background = element_blank(),axis.line = element_line(colour = "black"),
                     legend.position = "none",axis.text.y = element_text(size = font_size),
                     axis.text.x = element_text(size = font_size),#the font size of axis
                     axis.title.x = element_text(size = font_size),#the font size of axis title
                     axis.title.y = element_text(size = font_size))+
      ggplot2::geom_hline(yintercept= 1.30103,#draw a line at p=0.05
                          lty=4,col="grey21",lwd=1)+
      ggplot2::labs(title="Pathway Enrichment Analysis")+
      ggplot2::coord_flip()+
      ggplot2::xlab('Pathway')
    save(pb,file = 'barplot.Rda')
    export::graph2ppt(x=pb,file='pathway.pptx',height=7,width=9)
    dev.off()
    setwd('..//')
  }

}

