#' @title PathwayEnrich
#' @description a function to do PathwayAnalysis
#' @author Xia Shen
#' \email{qq951633542@@163.com}
#' @param font_size 20
#' @param specias_pathway_database the specia of sample
#' @param specias_compound_database the specia of sample
#' @return  All the results can be got form other functions and instruction.
#' @export
PathwayEnrich <- function(specias_pathway_database= c(hsa.kegg.pathway,mmu.kegg.pathway),
                          specias_compound_database = c(hsa_compound_ID,mmu_compound_ID),
                          font_size=20){
  pacman::p_load(ggplotify,readr,data.table,export,magrittr)
  data <- read_csv('data_pathway.csv')
  if(!file.exists('PathwayEnrich')){
    dir.create('PathwayEnrich')
  }

  #pathwayenrich
  #seed <- dplyr::filter(data,Annotation.type=='seed')
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
  # aID <- as.data.frame(IDinPathway[!duplicated(IDinPathway)])
  # i <- intersect(aID$`IDinPathway[!duplicated(IDinPathway)]`,specias_compound_database$id)
  # dt <- specias_compound_database[match(i,specias_compound_database$id),]
  # sd <- setdiff(aID$`IDinPathway[!duplicated(IDinPathway)]`,dt$id)
  # if(length(sd) != 0){
  #   nsd <- which(aID$`IDinPathway[!duplicated(IDinPathway)]`%in%sd)
  #   aa <- cbind(aID[-nsd,],dt)
  # }else(aa <- cbind(aID,dt))

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
  setwd('PathwayEnrich')
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
    scale_fill_manual(values = c('Turquoise3','Firebrick1'))+
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
  #####
  # vol<-read.csv("pathway_results.csv") #%>%
  #
  # Impact<- vol$Impact
  # p<- vol$Raw.p
  #
  #
  # sqx <- sqrt(Impact);
  # min.x<- min(sqx, na.rm = TRUE);
  # max.x <- max(sqx, na.rm = TRUE);
  #
  # if(min.x == max.x){ # only 1 value
  #   max.x = 1.5*max.x;
  #   min.x = 0.5*min.x;
  # }
  #
  # maxR <- (max.x - min.x)/40;
  # minR <- (max.x - min.x)/160;
  # radi.vec <- minR+(maxR-minR)*(sqx-min.x)/(max.x-min.x)
  # # set background color according to y
  # bg.vec <- heat.colors(length(p));
  # op <- par(mar=c(6,5,2,3));
  # plot(Impact, -log10(p), type="n", axes=F, xlab="Pathway Impact", ylab="-log10(p)");
  # axis(1);
  # axis(2);
  # grid(col="blue");
  # symbols(Impact, -log10(p), add = TRUE, inches = F, circles = radi.vec, bg = bg.vec, xpd=T);
  #
  # export::graph2ppt(file='pathway.pptx',height=7,width=9,append=T)

  ##
  # ph <- ggplot2::ggplot(vol, aes(x = Impact, y = -log10(Raw.p)))+
  #   geom_point(size=radi.vec*300,color=bg.vec) +
  #   #scale_color_manual(values = colour)
  #   # geom_vline(xintercept=0.75,
  #   #            lty=4,col="black",lwd=0.5)+ # add vetical line
  #   geom_hline(yintercept = -log10(0.05),
  #              lty=4,col="black",lwd=0.5)+ #add hori line
  #   labs(x="Pathway Impact",
  #        y="-log10(p-value)")+
  #   ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
  #                  panel.background = element_blank(),axis.line = element_line(colour = "black"),
  #                  legend.position = "none",axis.text.y = element_text(size = font_size),
  #                  axis.text.x = element_text(size = font_size),#the font size of axis
  #                  axis.title.x = element_text(size = font_size),#the font size of axis title
  #                  axis.title.y = element_text(size = font_size))+
  #   ggrepel::geom_text_repel(
  #     data = subset(vol,Raw.p< 0.05),
  #     #max.iter = 10,
  #     aes(label = X),
  #     size = path_font_size,
  #     box.padding = 0.25,
  #     point.padding = 0.3)
  # save(ph,file = 'pathway.Rda')
  # export::graph2ppt(x=ph,file='pathway.pptx',height=7,width=9,append=T)
  setwd('..//')
}

