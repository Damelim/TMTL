tab = matrix(NA, 2, 8)
colnames(tab) = c(paste('alpha = ', c(0,1/5,1/2,1)), paste('twosource', c('balanced','aligned')), paste('threesource', c('balanced','aligned')))
rownames(tab) = c('mean','sd')
library(dplyr)

load("alpha0/MTL_transfer/summary.Rdata")
if(length(iter_vec) == 20){
  tab[,1] = c(iter_vec %>% mean, iter_vec %>% sd)
}

load("alpha1:5/MTL_transfer/summary.Rdata")
if(length(iter_vec) == 20){
  tab[,2] = c(iter_vec %>% mean, iter_vec %>% sd)
}

load("alpha1:2/MTL_transfer/summary.Rdata")
if(length(iter_vec) == 20){
  tab[,3] = c(iter_vec %>% mean, iter_vec %>% sd)
}

load("alpha1/MTL_transfer/summary.Rdata")
if(length(iter_vec) == 20){
  tab[,4] = c(iter_vec %>% mean, iter_vec %>% sd)
}




load("twosource_balanced/MTL_transfer/summary.Rdata")
if(length(iter_vec) == 20){
  tab[,5] = c(iter_vec %>% mean, iter_vec %>% sd)
}

load("twosource_aligned/MTL_transfer/summary.Rdata")
if(length(iter_vec) == 20){
  tab[,6] = c(iter_vec %>% mean, iter_vec %>% sd)
}




load("threesource_balanced/MTL_transfer/summary.Rdata")
if(length(iter_vec) == 20){
  tab[,7] = c(iter_vec %>% mean, iter_vec %>% sd)
}

load("threesource_aligned/MTL_transfer/summary.Rdata")
if(length(iter_vec) == 20){
  tab[,8] = c(iter_vec %>% mean, iter_vec %>% sd)
}


tab %>% round(3)

