#' @title Calculate_ppm
#' @description a function to Calculate_ppm
#' @author Shine Shen
#' \email{qq951633542@@163.com}
#' @return  All the results can be got form other functions and instruction.
#' @export
Calculate_ppm <- function(){
  dat <- read_csv('rm_adducts.csv')
  a <- NULL
  for (i in seq(nrow(dat))) {
    b <- NULL
    l <- length(strsplit(dat$Formula,";")[[i]])
    bb <- c(b,l)
    a <- c(a,bb)
  }
  #ch3coo
  CH3COO <- dat[c(grep("CH3COO",dat$adduct),
                  which(dat$adduct == "")),]
  formula_ch3 <- separate(CH3COO[,'Formula'],Formula,into = letters[1:max(a)],sep = ";")
  f_ch3 <- apply(formula_ch3[,'a'], 1, function(x){
    ax <- paste(x,'CH3COO',sep = '')
  }) %>% as.data.frame()
  tmz_ch3 <- apply(f_ch3, 1, function(x){
    MonoisotopicMass(ListFormula(x),charge = 0)
  })
  ppm_ch3 <- abs((CH3COO$mz-tmz_ch3)/tmz_ch3*10^6 )
  mzc <- which(colnames(dat)=='mz')
  da_ch3 <- add_column(CH3COO,ppm_ch3,.after = 'mz')
  colnames(da_ch3)[mzc+1] <- 'ppm'

  #NH4
  NH4 <- dat[c(grep("NH4",dat$adduct),
               which(dat$adduct == "")),]
  formula_nh4 <- separate(NH4[,'Formula'],Formula,into = letters[1:max(a)],sep = ";")
  f_nh4 <- apply(formula_nh4[,'a'], 1, function(x){
    ax <- paste(x,'NH4',sep = '')
  }) %>% as.data.frame()
  tmz_nh4 <- apply(f_nh4, 1, function(x){
    MonoisotopicMass(ListFormula(x),charge = 0)
  })
  ppm_nh4 <- abs((NH4$mz-tmz_nh4)/tmz_nh4*10^6 )
  da_nh4 <- add_column(NH4,ppm_nh4,.after = 'mz')
  colnames(da_nh4)[mzc+1] <- 'ppm'
  #M-H
  MmH <- dat[c(grep("-H",dat$adduct),
               which(dat$adduct == "")),]
  formula_mh <- separate(MmH[,'Formula'],Formula,into = letters[1:max(a)],sep = ";")
  tmz_mh <- apply(formula_mh[,'a'], 1, function(x){
    MonoisotopicMass(ListFormula(x),charge = -1)
  })
  ppm_mh <- abs((MmH$mz-tmz_mh)/tmz_mh*10^6 )
  da_mh <- add_column(MmH,ppm_mh,.after = 'mz')
  colnames(da_mh)[mzc+1] <- 'ppm'
  #M+H
  MpH <- dat[-c(c(grep("CH3COO",dat$adduct),
                  which(dat$adduct == "")),c(grep("NH4",dat$adduct),
                                             which(dat$adduct == "")),c(grep("-H",dat$adduct),
                                                                        which(dat$adduct == ""))),]
  formula_ph <- separate(MpH[,'Formula'],Formula,into = letters[1:max(a)],sep = ";")

  tmz_ph <- apply(formula_ph, 1, function(x){
    MonoisotopicMass(ListFormula(x),charge = 1)
  })
  #
  ppm_ph <- abs((MpH$mz-tmz_ph)/tmz_ph*10^6 )
  da_ph <- add_column(MpH,ppm_ph,.after = 'mz')
  colnames(da_ph)[mzc+1] <- 'ppm'
  all <- rbind(da_ch3,da_nh4,da_mh,da_ph)
  write.csv(all,'ppm.csv')
}







