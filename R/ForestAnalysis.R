#' @title ForestAnalysis
#' @description a function can do forest analysis
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param group group set.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
ForestAnalysis<-function(group = c("case","control")){
  require(data.table)
  cat("Import data...\n")
  data <- fread("data.csv")
  data <- setDF(data)
  sample.info <- read.csv("sample.info.csv")

  case.name<-sample.info$sample.name[sample.info$group==group[1]]
  control.name<-sample.info$sample.name[sample.info$group==group[2]]

  case<-data[,match(case.name,colnames(data))]
  control<-data[,match(control.name,colnames(data))]

  mse <- apply(control,1,function(x) {
    u<-mean(x)+sd(x)
    l<-mean(x)-sd(x)
    c<-c(l,u)
  })

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
  require(meta)
  metaresult<-metabin(expose_1,case_all,expose_2,control_all,data=fdt,sm="OR",
                      studlab=paste(data$compound.name),comb.random=FALSE)
  tiff(file="forest plot.tiff", width = 1200, height = 1000,res = 56*2)
  forest(metaresult)
  dev.off()

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

  }))
  oddr <- data.frame(oddr[,1])
  colnames(oddr) <- "OR"

  dta <- cbind(dt,oddr)

  lowerci<-data.frame(apply(dta, 1, function(x){
    lower <- exp(log(dta[,5])-sqrt(1/dt[,1]+1/dt[,2]+1/dt[,3]+1/dt[,4])*1.96)

  }))
  lowerci <- data.frame(lowerci[,1])

  upperci<-data.frame(apply(dta, 1, function(x){
    upper <- exp(log(dta[,5])+sqrt(1/dt[,1]+1/dt[,2]+1/dt[,3]+1/dt[,4])*1.96)

  }))
  upperci <- data.frame(upperci[,1])
  result <- cbind(dta,lowerci,upperci)
  all <- cbind(oddr,lowerci,upperci,data)

  write.csv(all,"data_result.csv",row.names = F)

}
