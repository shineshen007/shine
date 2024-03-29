#' @title FeatureAnalysis
#' @description a function can generate mz vs rt plot and QC distribution plot,
#' also can filter isotope,rsd and zero value.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param zero.filter default is FALSE,if the zero value exsit in your data,
#' make FALSE to TRUE.
#' @param RSD.filter default is FALSE,if the percentage of qc rsd larger than 0.3
#'  more than 0.7,almost after normalization,make FALSE to TRUE.
#' @param zero.check default is TRUE.
#' @param mzrt default is FALSE.
#' @param FilterIsotope default is TRUE.
#' @return  All the results can be got form other functions and instruction.
#' @export
FeatureAnalysis <- function(zero.filter = FALSE,RSD.filter = FALSE,
                            zero.check = TRUE,mzrt = FALSE,FilterIsotope=TRUE) {

  cat("Import data...\n")
  data <- data.table::fread("data.csv")
  data <- data.table::setDF(data)
  sample.info <-read.csv("sample.info.csv")

  cat("Analyzing data...\n")
  ##create a folder for analysis
  dir.create("FeatureAnalysis")
  setwd("FeatureAnalysis")
  if(FilterIsotope){
  cat("Isotope filtering...\n")
  ###remove [M+n],\\ make [] lose the ability of function，
  isotope_filter<-function(data){
    temp<- data[c(grep("\\[M\\]",data$isotope),
                  which(data$isotope == "")),]
  }
  data <-isotope_filter(data)
  write.csv(data,"filter.isotope.csv",row.names = FALSE)
  }
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]#get sample index
  qc.name<-sample.info$sample.name[sample.info$class=="QC"]#get qc index

  sample<-data[,match(sample.name,colnames(data))]#get sample matrix
  qc<-data[,match(qc.name,colnames(data))]#get qc matrix

  sample.qc<-cbind(sample,qc)

  if(zero.check){
    cat("Zero checking...\n")
    ###zero check
    zero_check <- function(data){
      check_zero <- sample.qc
      numb.zero <- sapply(seq(nrow(check_zero)), function(i){
        check.num.zero <- sum(check_zero[i,]== 0)#get zero value columns
      })
      numb.sample <- ncol(check_zero)
      idx.check <- which(numb.zero >= numb.sample/2)
      check_zero <- check_zero[idx.check,]
    }
    zero.data <- zero_check(data)
    write.csv(zero.data,"zero.rows.csv",row.names = FALSE)
  }

  if(zero.filter){
    cat("Zero filtering...\n")
    ###zero filter
    zero_filter <- function(data){
      temp_zero <- data[,-c(1:4)]
      num.zero <- sapply(seq(nrow(temp_zero)), function(i){
        temp.num.zero <- sum(temp_zero[i,]== 0)
      })
      num.sample <- ncol(temp_zero)
      idx.filter <- which(num.zero >= num.sample/2)
      temp_zero <- temp_zero[-idx.filter,]
      temp <- data.frame(data[-idx.filter,c(1:4)], temp_zero)
    }
    data_zero_filter <- zero_filter(data)
    data <- data_zero_filter
    write.csv(data,"filter.zero.csv",row.names = FALSE)

  }

  cat("Calculate RSD...\n")
  ###calculate RSD
  RSD <- function(qc){
    qc_rsd <- qc
    rsd_qc <- sapply(seq(nrow(qc_rsd)), function(i){
      SD <- sd(qc_rsd[i,],na.rm = TRUE)
      MEAN <- sum(qc_rsd[i,],na.rm = TRUE)/ncol(qc)
      rsd_qc <- SD/MEAN
    })
  }
  rsd.data <- RSD(qc)##run function
  write.csv(rsd.data,"rsd.csv",row.names = FALSE)

  cat("Draw QC distribution plot...\n")
  ### QC distribution plot
  #png(file="QC distribution.png", width = 900, height = 800,res = 56*2)
  d <- read.csv("rsd.csv")
  percent <- round(sum(d$x<0.3)/nrow(d),3)
  txt <- paste(percent*100,"%")
  scatter.data<-as.data.frame(rsd.data[order(rsd.data)])
  id<-c(1:nrow(qc))
  data_qc<- cbind(scatter.data,id)
  qc_dis<- ggplot2::ggplot(data_qc,ggplot2::aes(x=id,y=rsd.data[order(rsd.data)]))+
    ggplot2::xlab("Feature Index")+
    ggplot2::ylab("Relative Standard Deviation(RSD)")+
    ggplot2::geom_point(ggplot2::aes(colour=scatter.data<0.3))+
    ggplot2::geom_text(data = d,aes(x= 400,y= 1.2,label= txt))
  ggplot2::geom_hline(ggplot2::aes(yintercept=0.3,linetype="dashed"))
  #plot(qc_dis)
  save(qc_dis,file = 'qc_dis.Rda')
  export::graph2ppt(x=qc_dis,file='qc_dis.pptx',height=7,width=9)
  #dev.off()

  if(RSD.filter){
    cat("RSD filtering...\n")
    ###RSD filter
    data_rsd <- data
    RSD_filter <- function(data_rsd){
      qc.rsd.name<-sample.info$sample.name[sample.info$class=="QC"]
      qc_rsd<-data_rsd[,match(qc.rsd.name,colnames(data_rsd))]
      temp_rsd <- qc_rsd
      rsd <- sapply(seq(nrow(temp_rsd)), function(i){
        SD <- sd(temp_rsd[i,],na.rm = TRUE)
        MEAN<-sum(temp_rsd[i,],na.rm = TRUE)/ncol(qc)
        rsd<-SD/MEAN
      })
      idx.filter <- which(rsd >= 0.3)
      rsd_data <- data_rsd[-idx.filter,]
    }
    filter.rsd.data <- RSD_filter(data_rsd)
    write.csv(filter.rsd.data,"data for sta.csv",row.names = FALSE)
  }
  cat("Draw mz VS RT plot...\n")
  if(mzrt){
    #### mz VS RT plot
    #png(file="mz.rt.png", width = 900, height = 800,res = 56*2)
    col <- apply(sample,1,median)
    mr <- ggplot2::ggplot(data,ggplot2::aes(x=rt,y=mz,colour=log10(col)))+
      geom_point()+
      scale_color_gradient(low = 'lightgreen', high = 'darkred')+
      ggplot2::xlab("Retention time(s)")+
      ggplot2::ylab("Mass to charge ratio(m/z)")+
      ggplot2::labs(colour="log10(intensity)")+
      ggplot2::theme(legend.position = c(0.95,0.9))
    #plot(mr)
    #dev.off()
    save(mr,file = 'mr.Rda')
    export::graph2ppt(x=mr,file='data.pptx',height=7,width=9)
  }

  ##back origin work directory
  setwd("..//")

  cat("FeatureAnalysis is done\n")
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage("Shine 0.3.7
                        Maintainer: Xia Shen.
                        \n2023-01-18
                        News: Update PathwayEnrich
                        Notes: sample name in pos and neg mode must be identical


                        Version 0.3.7
                        --------------
                        "
  )
}
