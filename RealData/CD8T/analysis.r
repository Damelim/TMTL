library(dplyr)
l2norm = function(x){
  sqrt(sum(x^2))
}

donorvec = 1:8 
n_donor_vec = rep(NA,length(donorvec))
load('clr_data_donor_nofiltered_1000hvg1.Rdata') ; K = ncol(Y) 
for(jj in 1:length(donorvec)){
  load(paste('clr_data_donor_nofiltered_1000hvg',jj,'.Rdata', sep = ""))
  n_donor_vec[jj] = nrow(Y)
}
n_donor_vec
rm(X,Y)
rm(donorvec, jj, K)


load(paste('MTL_transfer_debiased/summary_logtransform_target',1,'.Rdata',sep="")) ; nrep = length(W_opt_debiased_list)
num_common_features_target = rep(NA, nrep)
num_common_features_debiased =rep(NA,nrep)
top_common_features_target = list()
top_common_features_debiased=list()
for(k in 1:nrep){
  cat('=======================\n')
  cat('k=',k,'\n')
  cat('=======================\n')
  Btarget_list = list()
  for(j in 1:8){
    load(paste('MTL_notransfer/summary_logtransform_target',j,'.Rdata',sep=""))
    Btarget_list[[j]] = parameter_list[[k]]
  }
  Bdebiased_list = list()
  for(j in 1:8){
    load(paste('MTL_transfer_debiased/summary_logtransform_target',j,'.Rdata',sep=""))
    Bdebiased_list[[j]] = W_opt_debiased_list[[k]] ; rownames(Bdebiased_list[[j]]) = rownames(Btarget_list[[1]])
  }
  mat = matrix(NA,8,2) ; rownames(mat) = paste('donor',1:8,sep="") ; colnames(mat) = c('target','debiased')
  
  for(j in 1:8){
    #cat('j = ',j,'\n')
    mat[j,1] = 100*mean(rowMeans(as.matrix(Btarget_list[[j]])) == 0) ; mat[j,2] = 100*mean(rowMeans(as.matrix(Bdebiased_list[[j]])) == 0)
  }
  mat %>% round(2) %>% print
  
  selected_features_target = list() ; selected_features_debiased = list() 
  top_features_target = list() ; top_features_debiased = list()
  for(j in 1:8){
    #cat('j = ',j,'\n')
    selected_features_target[[j]] = names(which(rowMeans(as.matrix(Btarget_list[[j]])) != 0)) ; selected_features_debiased[[j]] = names(which(rowMeans(as.matrix(Bdebiased_list[[j]])) != 0))
    top_features_target[[j]] = names(sort(round(apply(Btarget_list[[j]], 1, l2norm),3), decreasing = T)[1:50])
    top_features_debiased[[j]] = names(sort(round(apply(Bdebiased_list[[j]], 1, l2norm),3), decreasing = T)[1:50])
  }
  Reduce(intersect, selected_features_target) %>% print ; 
  num_common_features_target[k] = Reduce(intersect, selected_features_target) %>% length
  Reduce(intersect, selected_features_debiased) %>% print 
  num_common_features_debiased[k] = Reduce(intersect, selected_features_debiased) %>% length
  
  top_common_features_target[[k]] = Reduce(intersect, top_features_target)
  top_common_features_debiased[[k]] = Reduce(intersect, top_features_debiased)
}
num_common_features_target %>% mean ; num_common_features_target %>% sd
num_common_features_debiased %>% mean ; num_common_features_debiased %>% sd

top_common_features_target
top_common_features_debiased
