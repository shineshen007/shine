#' @title Permutation_Shine
#' @description a function to do permutation test.
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  get data you want.
#' @param repeats the repeats to do,default is 200
#' @param ncomp the components of pls
#' @param group the group of the data
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
Permutation_Shine <- function(repeats = 200, ncomp = 2,
                           group = c('case','control')
                           ) {

  cat("Import data...\n")
  data <- data.table::fread("data.csv")
  data <- data.table::setDF(data)
  sample.info <- read.csv("sample.info.csv")
  class<- sample.info[,"group"]

  group1.index <- which(class == group[1])
  group2.index <- which(class == group[2])
  sample.info<-sample.info[c(group1.index,group2.index),]
  ###data preparation
  sample.name<-sample.info$sample.name[sample.info$class=="Subject"]
  sample<-data[,match(sample.name,colnames(data))]

  sample.index <- which(sample.info$class=="Subject")
  datat<-sample
  datatm<-as.matrix(datat)
  XXt<-t(datatm)
  group_pls<-as.data.frame(sample.info$group)
  YY<-as.numeric(group_pls[sample.index,])
  pls <- plsdepot::plsreg1(XXt, YY, comps = ncomp)

  Q2 <- pls$Q2[ncomp,5]
  R2 <- sum(pls$R2)

  ##begin repeat
  q2 <- NULL
  r2 <- NULL
  cor <- NULL
  for (i in 1:repeats) {
    Yt <- sample(YY)
    temp.pls <- plsdepot::plsreg1(XXt,Yt,comps = ncomp)
    q2[i] <- temp.pls$Q2[ncomp,5]
    r2[i] <- sum(temp.pls$R2)
    cor[i] <- abs(cor(Yt, YY))
    cat(i); cat(" ")
  }

  ##draw perumtation test
  FileName <- paste('permutation test', group[1], sep = " ")
  FileName <- paste(FileName, group[2], sep = "")
  FileName <- paste(FileName, ".pdf", sep = "")
  pdf(FileName)
  par(xpd = F)
  par(mar=c(5,5,4,2))
  plot(x = 0, y = 0, xlim = c(0,1),
       ylim = c(min(c(q2,r2,Q2,R2)),1.2*max(c(q2,r2,Q2,R2))),
       col = "white",
       xlab = "Correlation",
       ylab = "Values (Q2, R2)",
       cex.axis = 1.5, cex.lab = 1.8)
  abline( h = 0, lty = 2)

  points(x = cor, y = q2, col = "palegreen", pch = 19)
  points(x = cor, y = r2, col = "royalblue", pch = 19)

  points(x = 1, y = Q2, col = "palegreen", pch = 19)
  points(x = 1, y = R2, col = "royalblue", pch = 19)

  lm.r2 <- lm(c(R2,r2)~c(1, cor))
  lm.q2 <- lm(c(Q2,q2)~c(1, cor))

  intercept.q2 <- lm.q2$coefficients[1]
  intercept.r2 <- lm.r2$coefficients[1]

  segments(x0 = 0, y0 = intercept.q2, x1 = 1, y1 = Q2, lty = 2, lwd = 2)
  segments(x0 = 0, y0 = intercept.r2, x1 = 1, y1 = R2, lty =2, lwd = 2)


  legend("bottomright", title = "Intercepts",
         legend = c(paste("Q2",round(intercept.q2,2), sep = ": "), paste("R2",round(intercept.r2,2),sep=": ")),
         col = c("palegreen", "royalblue"), pch = 19, pt.cex = 1.3, cex = 1.3, bty = "n")
  par(xpd = T)
  dev.off()
}
