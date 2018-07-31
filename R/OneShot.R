#' @title OneShot
#' @description a function to do all work
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
OneShot<-function(path1=NULL,#dna/pos
                   path2=NULL,#dna/neg
                   path3=NULL,#dna both
                  group = c("case","control"),
                  pathway.enrichment= TRUE,
                  correct=TRUE
                   ){
  #First click: put peaktable.csv and sample.info.csv in pos and neg folder respectively
  #set working directory pos folder manually
  library(Shine)
  ArrangeData()#generate data for svr and info
  dir.create("svr")#in pos folder
  file.copy(c("data for svr.csv", "info.csv"), "svr")
  setwd("svr")#pos folder manually
  library(MetCleaning)
  MetCleaning(#ImportData para
    data = "data for svr.csv",
    sample.information = "info.csv",
    polarity = "positive",
    #DataNormalization
    method = "svr")
  DataID()#generate data for dna and info
  file.copy(c("data for dna.csv", "sample.info.csv"), path1)#need change

  setwd("..//")
  setwd("..//")
  setwd("neg")
  ArrangeData()#generate data for svr and info
  dir.create("svr")
  file.copy(c("data for svr.csv", "info.csv"), "svr")
  setwd("svr")
  library(MetCleaning)
  MetCleaning(#ImportData para
    data = "data for svr.csv",
    sample.information = "info.csv",
    polarity = "negative",
    #DataNormalization
    method = "svr")
  DataID()#generate data for svr and info
  file.copy(c("data for dna.csv", "sample.info.csv"), path2)#need change
  setwd("..//")
  setwd("..//")
  setwd("dna")
  library(MetDNA)#sample name in pos and neg must be identical
  ##正负合并处理
  MetDNA(ms1.data.pos = "data for dna.csv",
         ms1.data.neg = "data for dna.csv",
         sample.info.pos = "sample.info.csv",
         sample.info.neg = "sample.info.csv",
         pos.path = path1,#need change
         neg.path = path2,#need change
         polarity = "both",
         column = "hilic",
         ce = "30",
         use.default.md = TRUE,
         group = c("G","M"),
         uni.test = "t",
         correct = correct,
         p.cutoff = 0.05,
         species = "hsa",
         dn.analysis = FALSE,
         pathway.enrichment = pathway.enrichment)
  #pos file operation
  dir.create("both")
  setwd("both")
  setwd("..//")
  setwd("pos")
  file.copy("data for dna.csv", path3)#need change
  setwd("..//")
  setwd("both")
  file.rename("data for dna.csv","data pos.csv")
  setwd("..//")
  setwd("pos/MRN_annotation_result")
  file.copy("MRN.annotation.result.csv", path3)#need change
  setwd("..//")
  setwd("..//")
  setwd("both")
  file.rename("MRN.annotation.result.csv","ID POS.csv")
  #neg file operation
  setwd("..//")
  setwd("neg")
  file.copy("data for dna.csv", path3)#need change
  setwd("..//")
  setwd("both")
  file.rename("data for dna.csv","data neg.csv")
  setwd("..//")
  setwd("neg/MRN_annotation_result")
  file.copy("MRN.annotation.result.csv", path3)#need change
  setwd("..//")
  setwd("..//")
  setwd("both")
  file.rename("MRN.annotation.result.csv","ID NEG.csv")
  #
  BindData()
  data<-fread("data for next.csv")
  data<-setDF(data)
  NAfilter(data = data)
  da<-fread("filter.NA.csv")
  da<-setDF(da)
  write.csv(da,"data.csv",row.names = F)
  setwd("..//")
  setwd("..//")
  file.copy("sample.info.csv", path3)#need change
  setwd("dna/both")
  FeatureAnalysis(data = data,sample.info = sample.info,RSD.filter = T,
                  zero.check = F,zero.filter = F)
  setwd("FeatureAnalysis")
  dat<-fread("data for sta.csv")
  dat<-setDF(dat)
  dat<-dat[,-c(4,6,7)]
  write.csv(dat,"data final 5.csv",row.names = F)
  file.copy("data for sta.csv", path3)#need change
  setwd("..//")
  StaAnalysis(data = data,sample.info = sample.info,group = c("G","M"),
              pcorrect = T,unitest = t.test,p.cutoff = 0)
}




