#' @title BiClass
#' @description a function can do classify analysis using logistic,svm and randomforest
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
BiClass <- function(){
  #logistic regression
  require(pROC)
  require(e1071)
  require(randomForest)
  require(data.table)
  data <- fread("data.csv")
  data<-setDF(data)
  d <- dim(data)
  split <- sample(d[1],d[1]*(2/3))
  train<-data[split,]
  test<-data[-split,]
  fit.log<-glm(group~.,data = train,family = binomial())
  step.fit<-step(fit.log)
  prob<-predict(step.fit,test,type="response")
  logit.pred <- factor(prob > .5, levels=c(FALSE, TRUE),
                       labels=c("after", "before"))
  logit.perf <- table(test$group, logit.pred,
                      dnn=c("Actual", "Predicted"))
  tiff(file="logistic.tiff", width = 1200, height = 1000,res = 56*2)
  #logistic
  roc<-plot.roc(test[,1],prob,print.auc=T, max.auc.polygon=TRUE,print.thres=F,col="black")
  dev.off()
  #SVM

  fit.svm<-svm(group~.,data = train)
  svm.pred <- predict(fit.svm,test)
  svm.perf <- table(na.omit(test)$group,svm.pred, dnn=c("Actual", "Predicted"))

  pred.svm<-as.numeric(svm.pred)
  pred.svm[which(pred.svm==2)] <- 0
  tiff(file="svm.tiff", width = 1200, height = 1000,res = 56*2)
  #svm
  roc.svm<-plot.roc(test[,1],pred.svm,print.auc=T, print.thres=F,col="red")
  dev.off()

  #random forest

  fit.rf<-randomForest(group~.,data = train,na.action=na.roughfix,importance=TRUE)

  imp <- importance(fit.rf,type = 2)
  rf.pred<-predict(fit.rf,test)
  rf.perf<-table(test$group,rf.pred,dnn=c("Actual", "Predicted"))

  pred.rf<-as.numeric(rf.pred)
  pred.svm[which(pred.rf==2)] <- 0
  tiff(file="rf.tiff", width = 1200, height = 1000,res = 56*2)
  #rf
  roc.rf<-plot.roc(test[,1],pred.rf,print.auc=T,
                   col="blue",print.thres=F,auc.polygon=F)

  dev.off()
}

