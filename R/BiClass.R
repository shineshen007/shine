#' @title BiClass
#' @description a function can do classify analysis using logistic,svm and randomforest
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @param times the loop times in sample selection,default is 1001,and the number must be odd number.
#' @param logistic do logistic regression
#' @param RandomForest do RandomForest
#' @param SVM do support vector machine
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
BiClass <- function(times = 1001,#must be odd
                    logistic = TRUE,
                    RandomForest = TRUE,
                    SVM = TRUE
                    ){

  data <- utils::read.csv("data.csv",stringsAsFactors = T)
  d <- dim(data)
  num <- list()#save the split list
  au_lg <- NULL#save logistic auc
  au_svm <- NULL#save svm auc
  au_rf <- NULL#save rf auc
  set.seed(1234)
  for (i in 1:times ){
    split <- list(sample(d[1],d[1]*(2/3)))#split data
    numo <- list()#save the list of data selection
    nn <- c(numo,split)
    num <- c(num,nn)#save the list of split
    train<-data[num[[i]],]
    test<-data[-num[[i]],]
    #logistic
    fit.log <- glm(group~.,data = train,family = binomial())
    step.fit<-step(fit.log)
    prob<-predict(step.fit,test,type="response")
    roc_lg <- pROC::roc(test$group,prob,ci=T)
    auc_lg <- roc_lg[["auc"]][1]
    a_lg <- NULL
    aa_lg <- c(a_lg,auc_lg)
    au_lg <- c(au_lg,aa_lg)#get 1000 times auc value

    #svm
    fit.svm<-e1071::svm(group~.,data = train, probability = TRUE)
    svm.pred <- predict(fit.svm,test, probability = TRUE)
    pred.svm <- attr (svm.pred, "probabilities")[, 1]
    roc_svm <- pROC::roc(test[,1],pred.svm,ci=T)
    auc_svm <- roc_svm[["auc"]][1]
    a_svm <- NULL
    aa_svm <- c(a_svm,auc_svm)
    au_svm <- c(au_svm,aa_svm)

    #random forest

    fit.rf <- randomForest::randomForest(group~.,data = train,importance=TRUE, probability = TRUE)
    imp <- randomForest::importance(fit.rf,type = 2)
    rf.pred<-predict(fit.rf,test, type="prob")

    roc_rf <- roc(test[,1],rf.pred[,1],ci=T)
    auc_rf <- roc_rf[["auc"]][1]
    a_rf <- NULL
    aa_rf <- c(a_rf,auc_rf)
    au_rf <- c(au_rf,aa_rf)

  }
  #logistic
  if(logistic){
  AUC_lg_med <- median(au_lg)
  num_auc_lg <- which(au_lg==AUC_lg_med)[1]#select the median auc split

  ind_lg <- num[[num_auc_lg]]#get the median auc split
  train_data_lg<-data[ind_lg,]
  test_data_lg<-data[-ind_lg,]
  fit.log<-glm(group~.,data = train_data_lg,family = binomial())
  step.fit<-step(fit.log)
  prob<-predict(step.fit,test_data_lg,type="response")
  roc_lg <- roc(test_data_lg[,1],prob, ci=T)

  #
  upper_lg <- round(roc_lg[["ci"]][3],2)
  med_lg <- round(median(au_lg),2)
  lower_lg <- round(roc_lg[["ci"]][1],2)

  cc_lg <- data.frame(rep("1",times))
  v <- cbind(cc_lg,au_lg)
  colnames(v) <- c("index","auc")
  label="95% CI"
  p <- ggplot2::ggplot(mapping= ggplot2::aes(v$index, v$auc)) +
    ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),#remove ggplot2 background
          panel.background = element_blank(),axis.line = ggplot2::element_line(colour = "black"),legend.position = "none")+
    ggbeeswarm::geom_quasirandom(aes(color="grey"))+
    labs(x=NULL,
         y="Area Under Curve(AUC)",
         title="AUC Distribution plot")+
    ggplot2::geom_segment(mapping = aes(x=0.9,y=upper_lg),xend=1.1,yend=upper_lg)+
    geom_segment(mapping = aes(x=0.7,y=med_lg),xend=1.3,yend=med_lg)+
    geom_segment(mapping = aes(x=0.9,y=lower_lg),xend=1.1,yend=lower_lg)+
    geom_segment(mapping = aes(x=1,y=upper_lg),xend=1,yend=lower_lg)+
    geom_segment(mapping = aes(x=1.4,y=upper_lg),xend=1.4,yend=lower_lg,lty=4)+
    annotate("text", x=1.45, y=med_lg, label=label, colour='black', size=4)+
    annotate("text", x=1.25, y=upper_lg-0.01, label=upper_lg, colour='black', size=5)+
    annotate("text", x=1.25, y=med_lg-0.01, label=med_lg, colour='black', size=5)+
    annotate("text", x=1.25, y=lower_lg-0.01, label=lower_lg, colour='black', size=5)+
    ggplot2::ggsave("AUC logistic.png",width=10,height=6)
  }

  #SVM
  if(SVM){
  AUC_svm_med <- median(au_svm)
  num_auc_svm <- which(au_svm==AUC_svm_med)[1]

  ind <- num[[num_auc_svm]]
  train_data_svm<-data[ind,]
  test_data_svm<-data[-ind,]
  fit.svm<-e1071::svm(group~.,data = train_data_svm,probability = T)
  svm.pred <- predict(fit.svm,test_data_svm, type="prob",probability = T)
  prob.svm <- attr (svm.pred, "probabilities")[, 1]
  roc_svm <- pROC::roc(test_data_svm[,1],prob.svm, ci=T)
  auc_svm <- roc_svm[["auc"]][1]
  upper_svm <- round(roc_svm[["ci"]][3],2)
  med_svm <- round(median(au_svm),2)
  lower_svm <- round(roc_svm[["ci"]][1],2)

  cc <- data.frame(rep("1",times))
  v <- cbind(cc,au_svm)
  colnames(v) <- c("index","auc")
  label="95% CI"
  p <- ggplot(mapping=aes(v$index, v$auc)) +
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.position = "none")+
    ggbeeswarm::geom_quasirandom(aes(color="grey"))+
    labs(x=NULL,
         y="Area Under Curve(AUC)",
         title="AUC Distribution plot")+
    geom_segment(mapping = aes(x=0.9,y=upper_svm),xend=1.1,yend=upper_svm)+
    geom_segment(mapping = aes(x=0.7,y=med_svm),xend=1.3,yend=med_svm)+
    geom_segment(mapping = aes(x=0.9,y=lower_svm),xend=1.1,yend=lower_svm)+
    geom_segment(mapping = aes(x=1,y=upper_svm),xend=1,yend=lower_svm)+
    geom_segment(mapping = aes(x=1.4,y=upper_svm),xend=1.4,yend=lower_svm,lty=4)+
    annotate("text", x=1.45, y=med_svm, label=label, colour='black', size=4)+
    annotate("text", x=1.25, y=upper_svm-0.01, label=upper_svm, colour='black', size=5)+
    annotate("text", x=1.25, y=med_svm-0.01, label=med_svm, colour='black', size=5)+
    annotate("text", x=1.25, y=lower_svm-0.01, label=lower_svm, colour='black', size=5)+
    ggsave("AUC svm.png",width=10,height=6)

}

  #random forest
  if(RandomForest){
  AUC_rf_med <- median(au_rf)
  num_auc_rf <- which(au_rf==AUC_rf_med)[1]
  ind <- num[[num_auc_rf]]
  train_data<-data[ind,]
  test_data<-data[-ind,]
  fit.rf<-randomForest(group~.,data = train_data,importance=TRUE, probability = TRUE)
  imp <- importance(fit.rf,type = 2)
  rf.pred<-predict(fit.rf,test_data, type="prob")

  roc_rf <- pROC::roc(test_data[,1],rf.pred[,1], ci=T)

  upper_rf <- round(roc_rf[["ci"]][3],2)
  med_rf <- round(median(au_rf),2)
  lower_rf <- round(roc_rf[["ci"]][1],2)
  png(file="ROC.png", width = 1200, height = 1000,res = 56*2)
  #roclg
  rocc_lg<-pROC::plot.roc(roc_lg,col="black"#, print.auc = T
                    #,print.thres = "best"
  )
  #rocsvm
  rocc<-plot.roc(roc_svm,col="blue",add = TRUE#,print.auc = T
                 #,print.thres = "best"
  )
  #rocrf
  rocc<-plot.roc(roc_rf,col="red",add = TRUE#,print.auc = T
                 #,print.thres = "best"
  )

  legend(0.65,0.15, legend=c('Logistic Regression', "SVM","RandomForest"),
         col=c("black", "blue","red"), lwd=2,bty = "n")
  legend(0.4,0.15, legend=c(med_lg,med_svm,med_rf),
         col=c("black", "blue","red"), lwd=0,bty = "n")
  legend(0.3,0.15, legend=c(upper_lg,upper_svm,upper_rf),
         col=c("black", "blue","red"), lwd=0,bty = "n")
  legend(0.24,0.15, legend=c('--',"--","--"),
         col=c("black", "blue","red"), lwd=0,bty = "n")
  legend(0.21,0.15, legend=c(lower_lg,lower_svm,lower_rf),
         col=c("black", "blue","red"), lwd=0,bty = "n")
  text(0.22,0.15,labels = "AUC (95%CI)")
  dev.off()

  cc <- data.frame(rep("1",times))
  v <- cbind(cc,au_rf)
  colnames(v) <- c("index","auc")
  label="95% CI"
  p <- ggplot2::ggplot(mapping = ggplot2::aes(v$index, v$auc)) +
    ggplot2::theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.position = "none")+
    geom_quasirandom(ggplot2::aes(color="grey"))+
    labs(x=NULL,
         y="Area Under Curve(AUC)",
         title="AUC Distribution plot")+
    geom_segment(mapping = aes(x=0.9,y=upper_rf),xend=1.1,yend=upper_rf)+
    geom_segment(mapping = aes(x=0.7,y=med_rf),xend=1.3,yend=med_rf)+
    geom_segment(mapping = aes(x=0.9,y=lower_rf),xend=1.1,yend=lower_rf)+
    geom_segment(mapping = aes(x=1,y=upper_rf),xend=1,yend=lower_rf)+
    geom_segment(mapping = aes(x=1.4,y=upper_rf),xend=1.4,yend=lower_rf,lty=4)+
    annotate("text", x=1.45, y=med_rf, label=label, colour='black', size=4)+
    annotate("text", x=1.25, y=upper_rf-0.01, label=upper_rf, colour='black', size=5)+
    annotate("text", x=1.25, y=med_rf-0.01, label=med_rf, colour='black', size=5)+
    annotate("text", x=1.25, y=lower_rf-0.01, label=lower_rf, colour='black', size=5)+
    ggsave("AUC rf.png",width=10,height=6)
}
}

