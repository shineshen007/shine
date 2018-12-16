#' @title BoxPlot
#' @description a function to draw BoxPlot
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @param metabolite the metabolite to build boxplot
#' @export
#' @examples
#' \donttest{
#' ##---- Be sure the format of data and sample.info is correct!! ----
#' }
BoxPlot <- function(metabolite=NULL){

  data <- read.csv("data.csv",check.names = F)
  #gender
  Gender <- as.data.frame(data$sex)

  colnames(Gender)<-"Gender"

  col <- colnames(data)
  m <- match(metabolite,col)
  d <- data[,m]
  scale <- scale(d)
  ylab <- metabolite
  data1 <- cbind(data,Gender,scale)
  tiff(file="boxplot of sex.tiff", width = 1200, height = 700,res = 56*2)
  s <- ggpubr::ggboxplot(data = data1,x="Gender",y = "scale",color = "Gender",
            palette = "jco",add = "jitter",ylab = ylab)+
  stat_compare_means()
  plot(s)
  dev.off()

  #diebite

  tiff(file="boxplot of diebite.tiff", width = 1200, height = 700,res = 56*2)
  diabetes <- ggboxplot(data = data1,x="diabetes",y = "scale",color = "diabetes",
                 palette = "jco",add = "jitter",ylab = ylab)+
    stat_compare_means()
  plot(diabetes)
  dev.off()

  #hypertension

  tiff(file="boxplot of hypertension.tiff", width = 1200, height = 700,res = 56*2)
  hypertension <- ggboxplot(data = data1,x="hypertension",y = "scale",color = "hypertension",
                 palette = "jco",add = "jitter",ylab = ylab)+
    stat_compare_means()
  plot(hypertension)
  dev.off()

  # hyperlipemia

  tiff(file="boxplot of hyperlipemia.tiff", width = 1200, height = 700,res = 56*2)
  hyperlipemia <- ggboxplot(data = data1,x="hyperlipemia",y = "scale",color = "hyperlipemia",
                 palette = "jco",add = "jitter",ylab = ylab)+
    stat_compare_means()
  plot(hyperlipemia)
  dev.off()

  #smoke
  tiff(file="boxplot of smoke.tiff", width = 1200, height = 700,res = 56*2)
  smoke <- ggboxplot(data = data1,x="smoke",y = "scale",color = "smoke",
                 palette = "jco",add = "jitter",ylab = ylab)+
    stat_compare_means()
  plot(smoke)
  dev.off()

  #age
  age <- table(cut(data$age, breaks = c(0,49,60,61,70,max(data$age))))
  age <- as.vector(age)
  a <- age[1]
  b <- age[2]
  c <- age[3]+age[4]
  d <- age[5]
  Age <- as.data.frame(c(rep("<50",a),rep("50-60",b),rep("61-70",c),rep("70+",d)))
  colnames(Age)<-"Age"
  data <- data[order(data[,"age"]),] #sort by age
  data2 <- cbind(data,Age,scale)
  my_comparisons <- list(c("<50", "50-60"), c("<50", "61-70"), c("<50", "70+"),c("50-60","61-70")
                         ,c("50-60","70+"),c("61-70","70+"))
  tiff(file="boxplot of age.tiff", width = 1200, height = 700,res = 56*2)
  age <- ggboxplot(data = data2,x="Age",y = "scale",color = "Age",xlab = "years",
            palette = "jco",add = "jitter",ylab = ylab)+
    stat_compare_means(comparisons=my_comparisons,label.y = c(29, 35, 40,45,50,55))+
    stat_compare_means(label.y = 45)
  plot(age)
  dev.off()
}

