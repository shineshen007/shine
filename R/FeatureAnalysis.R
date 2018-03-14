#' @title FeatureAnalysis
#' @description a function can generate mz vs rt plot and QC distribution plot,
#' also can filter isotope,rsd and zero value.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param data a dataframe include name,mz,rt and isotope columns,
#' the rest of all are sample and QC columns..
#' @param sample.info a dataframe include sample.name,injection.order,
#' class,batch and group columns.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#load the demo data
#'data(data, package = "Shine")
#'data(sample.info, package = "Shine")

#'##create a folder for demo
#'dir.create("demo")
#'setwd("demo")

#'# export the demo data as csv
#'write.csv(data, "data.csv",sep = ",", row.names = FALSE)
#'write.csv(sample.info, "sample.info.csv",sep = ",", row.names = FALSE)

#'# Analysis process
#'FeatureAnalysis(data = "data.csv",sample.info = "sample.info.csv")
#' }
FeatureAnalysis <- function(data = NULL,sample.info = NULL) {
  cat("Analyzing data...\n")
  cat("Isotope filtering...\n")
###remove [M+n],\\ make [] lose the ability of function，
isotop_filter<-function(data){
  temp<- data[c(grep("\\[M\\]",data$isotope),
                which(data$isotope == "")),]
}
data <-isotop_filter(data)
write.csv(data,"filter.isotope.csv",row.names = FALSE)

###data preparation

  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  qc.name<-sample.info$sample.name[sample.info$class=="QC"]

  sample<-data[,match(sample.name,colnames(data))]
  qc<-data[,match(qc.name,colnames(data))]

  tags <- data[,c(1:4)]

  sample.tag<-cbind(tags,sample)

  rownames(sample)<-rownames(qc)<-tags$name

  cat("Calculate RSD...\n")
###calculate RSD
  RSD <- function(qc){
    temp_rsd <- qc
    rsd <- sapply(seq(nrow(temp_rsd)), function(i){
      SD <- sd(temp_rsd[i,])
      MEAN<-sum(temp_rsd[i,])/ncol(qc)
      rsd<-SD/MEAN
    })
  }
  rsd.data <- RSD(qc)##run function
  write.csv(rsd.data,"rsd.csv",row.names = FALSE)

  cat("Draw QC distribution plot...\n")
### QC distribution plot
  png(file="QC distribution.png", width = 900, height = 800,res = 56*2)
  scatter.data<-as.data.frame(rsd.data[order(rsd.data)])
  id<-c(1:nrow(data))
  data_qc<- cbind(scatter.data,id)
  qc_dis<- ggplot(data_qc,aes(x=id,y=rsd.data[order(rsd.data)]))+
    xlab("Feature Index")+
    ylab("Relative Standard Deviation(RSD)")+
    geom_point(aes(colour=scatter.data<0.3))+
    geom_hline(aes(yintercept=0.3,linetype="dashed"))
  plot(qc_dis)
  dev.off()

  cat("RSD filtering...\n")
  ###RSD filter
  data_rsd<-cbind(sample.tag,qc)
  RSD_filter <- function(data_rsd){
    temp_rsd <- data_rsd[,-c(1:ncol(sample.tag))]
    rsd <- sapply(seq(nrow(temp_rsd)), function(i){
      SD <- sd(temp_rsd[i,])
      MEAN<-sum(temp_rsd[i,])/ncol(qc)
      rsd<-SD/MEAN})
   idx.filter <- which(rsd >= 0.3)
   temp_rsd <- temp_rsd[-idx.filter,]
   temp_rsd <- data.frame(data_rsd[-idx.filter,c(1:ncol(sample.tag))], temp_rsd)
    }
  filter.rsd.data <- RSD_filter(data_rsd)
  write.csv(filter.rsd.data,"filter.rsd.csv",row.names = FALSE)

  cat("Zero filtering...\n")
  ###zero filter
  zero_filter <- function(filter.rsd.data){
    temp_zero <- filter.rsd.data
    num.zero <- sapply(seq(nrow(temp_zero)), function(i){
    temp.num.zero <- sum(temp_zero[i,]== 0)
  })
  num.sample <- ncol(temp_zero)
  idx.filter <- which(num.zero >= num.sample/2)
  temp_zero <- temp_zero[-idx.filter,]
  temp_zero <- data.frame(filter.rsd.data[-idx.filter,], temp_zero)
}
filter.zero.data <- zero_filter(filter.rsd.data)
write.csv(filter.zero.data,"filter.zero.csv",row.names = FALSE)

cat("Draw mz VS RT plot...\n")
#### mz VS RT plot
png(file="mz.rt.png", width = 900, height = 800,res = 56*2)
col <- apply(data[,-c(1:4,ncol(qc))],1,median)
mr <- ggplot(data,aes(x=rt,y=mz,colour=log10(col)))+
  geom_point()+
  scale_color_gradient(low = 'lightgreen', high = 'darkred')+
  xlab("Retention time(s)")+
  ylab("Mass to charge ratio(m/z)")+
  labs(colour="log10(intensity)")+
  theme(legend.position = c(0.95,0.9))
plot(mr)
dev.off()

cat("FeatureAnalysis is done\n")
}
