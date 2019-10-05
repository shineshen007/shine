#' @title CoxAnalysis
#' @description a function to do cox Proportional-Hazards Model
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
CoxAnalysis <- function(){
  data<- read.csv("data.csv")#data contain columns of patients'name,age,sex etc

  ###unvariate analysis
  cln <- colnames(data)
  len <- length(cln)
  covariates<-cln[4:len]
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(time, status)~', x)))

  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})
  # Extract data
  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (",
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  as.data.frame(res)
  write.csv(res,"uniresult.csv",row.names = T)

  ###multivariate cox
  covariatesp <- paste0(covariates[1:9], " +")
  clen <- length(covariates)
  x <- c(covariatesp,covariates[clen])
  formulas <- as.formula(paste('Surv(time, status)~', paste(x, collapse= "+")))
  mres.cox <- survival::coxph(formulas, data = data)
  ##forest plot
  #png(file="forest plot of multiv_cox.png", width = 1200, height = 1000,res = 56*2)
  sf <- survminer::ggforest(mres.cox,main = "Hazard ratio",cpositions = c(0.02, 0.22, 0.4),
           fontsize = 0.8,refLabel = "reference", noDigits = 2,data = data)
  export::graph2ppt(x=sf,file='correlation.pptx',height=7,width=9)
  #dev.off()
  summary(mres.cox)

  x <- summary(mres.cox)
  p.value<-signif(x$coefficients[,5], digits=2)
  HR <-signif(x$conf.int[,1], digits=2);#exp(beta)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
  HR <- paste0(HR, " (",
               HR.confint.lower, "-", HR.confint.upper, ")")
  HR<-as.data.frame(HR)
  p.value<-as.data.frame(p.value)
  res<-cbind(HR, p.value)
  names(res)<-c("HR (95% CI for HR)",
                "p.value")
  write.csv(res,"multi_result.csv",row.names = T)
}

