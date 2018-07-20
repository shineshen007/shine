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
png(file="forest plot.png", width = 1200, height = 1000,res = 56*2)
forest(metaresult)
dev.off()
}
