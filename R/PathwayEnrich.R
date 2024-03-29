#' @title PathwayEnrich
#' @description a function to do PathwayAnalysis
#' @author Xia Shen
#' \email{qq951633542@@163.com}
#' @param font_size 20
#' @param seed default is TRUE
#' @param All_ID default is FALSE
#' @param grade1_2 default is FALSE
#' @param data_name the data_name
#' @param specias_pathway_database the specia of sample
#' @param specias_compound_database the specia of sample
#' @param color color
#' @return  All the results can be got form other functions and instruction.
#' @export
PathwayEnrich <- function(specias_pathway_database= hsa.kegg.pathway,
                          specias_compound_database = hsa_compound_ID,
                          font_size = 20,
                          data_name = 'data_pathway.csv',
                          color = c("#3399FF","#FF3333"),
                          All_ID = TRUE,
                          seed = FALSE,
                          grade1_2 = FALSE
){
  pacman::p_load(ggplotify,readr,crayon,export,magrittr)
  if(All_ID){
    data <- read_csv(data_name) %>% dplyr::filter(fc>1)
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
    ##
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

    info <- info[order(info[,1]),]
    write.csv(info,'PathwayEnrich_All_ID/pathway_results.csv')
    #o('pathway_impact_76_pathway_hsa.all.ID.csv')
    sig_path <- dplyr::filter(info,q.value<0.05)
    if(nrow(sig_path)==0){
      cat(red('No significant pathway in enrichment analysis.\n'))
    }else{
      ip <- lapply(specias_pathway_database, function(x){
        intersect(x,mid)
      })

      for (i in setdiff(names(specias_pathway_database),rownames(sig_path))) {
        ip[[i]]=NULL
      }

      for (i in 1:length(ip)) {
        ip[[i]] <- data[match(ip[[i]],data$ID),]
      }
      sheetnames <- nchar(names(ip))
      if(max(sheetnames)>31){
        names(ip)[which(sheetnames>31)]=substr(names(ip)[which(sheetnames>31)],1,30)
      }
      setwd('PathwayEnrich_All_ID')
      openxlsx::write.xlsx(ip,'Id_pathway.xlsx',overwrite = T)
      #
      data <- read.csv("pathway_results.csv")
      colnames(data)[3] <- 'p'
      row=sum(data$p<0.05)+3
      data <- data[1:row,]
      group <- ifelse(data$p < 0.05,"sig", "not sig")
      pb <- ggplot2::ggplot(data,ggplot2::aes(reorder(X,-p),-log10(p)))+##-p control the order
        ggplot2::geom_bar(aes(fill=group),stat = "identity",position="dodge",width=0.8)+
        scale_fill_manual(values = color)+
        ggplot2::theme(panel.grid.major =element_blank(),
                       panel.grid.minor = element_blank(),#remove ggplot2 background
                       panel.background = element_blank(),
                       plot.title = element_text(size=24),
                       axis.line = element_line(colour = "black"),
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
      ggsave('pathway_barplot.pdf',width = 12,height = 9)
      #export::graph2ppt(x=pb,file='pathway.pptx',height=7,width=9)
      ####colnames(data)
      colnames(data) <- c('pathway','p.value','p.adjust','FDR','Pathway.length','Count')
      #####need modify
      d1=data %>% dplyr::filter(p.adjust<0.05&Count>1)
      d1=d1[order(d1$p.adjust),]
      d1$pathway <- fct_inorder(d1$pathway)
      d1$Pathway_impact=d1$Count/d1$Pathway.length
      #d1$Count=d1$Count*2
      dpp <- ggplot2::ggplot(d1, aes(x = Pathway_impact,y = pathway,
                                     color=p.adjust,size=Count))+
        geom_point()+
        #scale_fill_brewer(palette = "Set1")
        scale_color_continuous(low="red", high="blue",
                               guide=guide_colorbar(reverse=TRUE)) +
        theme_bw()+
        ggplot2::theme(#panel.grid.major =element_blank(),
          #                panel.grid.minor = element_blank(),#remove ggplot2 background
          #                panel.background = element_blank(),
          plot.title = element_text(size=24),
          legend.position = "right",
          axis.text.y = element_text(size = font_size),
          axis.text.x = element_text(size = font_size),#the font size of axis
          axis.title.x = element_text(size = font_size),#the font size of axis title
          axis.title.y = element_text(size = 0))
      ggsave('enrichKEGG.pdf',width = 12,height = 9)
      save(dpp,file = 'enrichKEGG_dotplot.RData')
      setwd('..//')
    }

  }
  #####
  if(seed){
    dir.create('PathwayEnrich_seed')
    #pathwayenrich
    data <- read_csv(data_name)
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

    info <- info[order(info[,1]),]
    write.csv(info,'PathwayEnrich_seed/pathway_results.csv')
    #o('pathway_impact_76_pathway_hsa.all.ID.csv')
    ##id match to pathway
    #compound_db <- read_csv('/Users/shineshen/Downloads/demo/kegg/old/kegg_compound_unique.csv')
    sig_path <- dplyr::filter(info,q.value<0.05)
    if(nrow(sig_path)==0){
      cat(red('No significant pathway in enrichment analysis.\n'))
    }else{
      ip <- lapply(specias_pathway_database, function(x){
        intersect(x,mid)
      })

      for (i in setdiff(names(specias_pathway_database),rownames(sig_path))) {
        ip[[i]]=NULL
      }

      for (i in 1:length(ip)) {
        ip[[i]] <- data[match(ip[[i]],data$ID),]
      }
      sheetnames <- nchar(names(ip))
      if(max(sheetnames)>31){
        names(ip)[which(sheetnames>31)]=substr(names(ip)[which(sheetnames>31)],1,30)
      }
      setwd('PathwayEnrich_seed')
      openxlsx::write.xlsx(ip,'Id_pathway.xlsx',overwrite = T)
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
      ggsave('pathway_barplot.pdf',width = 12,height = 9)
      #dev.off()
      setwd('..//')
    }}
  ###########
  if(grade1_2){
    dir.create('PathwayEnrich_grade1_2')
    #pathwayenrich
    #data <- dplyr::filter(data,Annotation.type=='seed')
    data <- read_csv(data_name)
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

    info <- info[order(info[,1]),]
    write.csv(info,'PathwayEnrich_grade1_2/pathway_results.csv')
    #o('pathway_impact_76_pathway_hsa.all.ID.csv')
    sig_path <- dplyr::filter(info,q.value<0.05)
    if(nrow(sig_path)==0){
      cat(red('No significant pathway in enrichment analysis.\n'))
    }else{
      ip <- lapply(specias_pathway_database, function(x){
        intersect(x,mid)
      })

      for (i in setdiff(names(specias_pathway_database),rownames(sig_path))) {
        ip[[i]]=NULL
      }

      for (i in 1:length(ip)) {
        ip[[i]] <- data[match(ip[[i]],data$ID),]
      }
      sheetnames <- nchar(names(ip))
      if(max(sheetnames)>31){
        names(ip)[which(sheetnames>31)]=substr(names(ip)[which(sheetnames>31)],1,30)
      }
      setwd('PathwayEnrich_grade1_2')
      openxlsx::write.xlsx(ip,'Id_pathway.xlsx',overwrite = T)
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
      ggsave('pathway_barplot.pdf',width = 12,height = 9)
      setwd('..//')
    }
  }
}

