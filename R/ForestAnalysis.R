#' @title ForestAnalysis_Shine
#' @description a function can do forest analysis
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param group group set.
#' @return  All the results can be got form other functions and instruction.
#' @import meta
#' @import survminer
#' @import survival
ForestAnalysis_Shine <-function(group = c("Gout","Normal")){

  cat("Import data...\n")
  data <- data.table::fread("data.csv")
  data <- data.table::setDF(data)
  sample.info <- read.csv("sample.info.csv")

  case.name<-sample.info$sample.name[sample.info$group==group[1]]#get case name
  control.name<-sample.info$sample.name[sample.info$group==group[2]]
  qc.name<-sample.info$sample.name[sample.info$class=="QC"]#get qc index
  qc.idx<-match(qc.name,colnames(data))#get qc index
  data <- data[,-qc.idx]

  case<-data[,match(case.name,colnames(data))]
  control<-data[,match(control.name,colnames(data))]

  mse <- apply(control,1,function(x) {
    u<-mean(x)+sd(x)
    l<-mean(x)-sd(x)
    c<-c(l,u)
  })#caculate the cutoff

  a<-apply(case, 1, function(x){
    length(which(x<mse[1]))
  })
  b<-apply(case, 1, function(x){
    length(which(x>mse[2]))
  })
  expose_1<-a+b
  case_all<-length(case)

  c<-apply(control, 1, function(x){
    length(which(x<mse[1]))
  })
  d<-apply(control, 1, function(x){
    length(which(x>mse[2]))
  })
  expose_2<-c+d
  control_all<-length(control)

  fd<-cbind(expose_1,case_all,expose_2,control_all)
  p<-as.data.frame(data$compound.name)
  fdt<-as.data.frame(cbind(p,fd))

  metaresult<- meta::metabin(expose_1,case_all,expose_2,control_all,data=fdt,sm="OR",
                             studlab=paste(data$name),comb.random=FALSE)
  #png(file="forest plot.png", width = 1200, height = 1000,res = 56*2)
  mf <- meta::forest(metaresult,leftlabs = c("Metabolites",NA,NA,NA,NA))
  save(mf,file = 'forest.Rda')
  export::graph2ppt(x=mf,file='forest plot.pptx',height=7,width=9)
  #dev.off()

  e <- as.data.frame(expose_1)
  f <- as.data.frame(expose_2)
  n <- length(expose_1)
  j<-data.frame(NULL)
  for (i in 1:n) {
    k<-NULL
    g <- as.data.frame(case_all-expose_1[i])
    temp.c <- rbind(k,g)
    j <- rbind(j, temp.c)
  }

  q<-data.frame(NULL)
  for (i in 1:n) {
    k<-NULL
    h <- as.data.frame(control_all-expose_2[i])
    temp.c <- rbind(k,h)
    q <- rbind(q, temp.c)
  }

  dt <- cbind(e,f,j,q)

  oddr<-data.frame(apply(dt, 1, function(x){
    or <- (dt[,1]/dt[,3])/(dt[,2]/dt[,4])

  }))#caculate OR
  oddr <- data.frame(oddr[,1])
  colnames(oddr) <- "OR"

  dta <- cbind(dt,oddr)

  lowerci<-data.frame(apply(dta, 1, function(x){
    lower <- exp(log(dta[,5])-sqrt(1/dt[,1]+1/dt[,2]+1/dt[,3]+1/dt[,4])*1.96)

  }))#caculate 95% CI
  lowerci <- data.frame(lowerci[,1])

  upperci<-data.frame(apply(dta, 1, function(x){
    upper <- exp(log(dta[,5])+sqrt(1/dt[,1]+1/dt[,2]+1/dt[,3]+1/dt[,4])*1.96)

  }))
  upperci <- data.frame(upperci[,1])
  result <- cbind(dta,lowerci,upperci)
  all <- cbind(oddr,lowerci,upperci,data)
  write.csv(all,"data_result.csv",row.names = F)

  filterOR <- all[!is.na(all$OR),]
  filterlow <- filterOR[!is.na(filterOR$lowerci...1.),]
  filterupper <- filterlow[!is.na(filterlow$upperci...1.),]
  fn <- which(filterupper$lowerci...1.<1&filterupper$upperci...1.>1)
  droc <- filterupper[-fn,]
  write.csv(droc,"data sig_OR.csv",row.names = F)
  #
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  rd<-droc[,match(sample.name,colnames(droc))]%>%
    t(.)
  colnames(rd) <- droc$name
  rownames(rd)=sample.info$group[match(rownames(rd),sample.info$sample.name)]
  #
  cs <- which(rownames(rd) == group[1])
  ct <- which(rownames(rd) == group[2])
  as <- rd[c(cs,ct),]
  write.csv(as,'ls.csv')
  roc <- read_csv('ls.csv')
  colnames(roc)[1] <- 'group'
  unlink('ls.csv')
  write.csv(roc,"data for roc.csv",row.names = F)
}
