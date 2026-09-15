library(glmnet)
library(Matrix)
library(dplyr)
source('functions_joint.R')
source('simulationsetting.R')

tgtidx = 5

load(paste('../clr_data_donor_nofiltered_500hvg',tgtidx,'.Rdata',sep=""))
X0 = X ; Y0 = Y ; Y0 = as.matrix(Y0) 
rm(X, Y) 

#### domain-wise, training-only centering of the predictors and the CLR responses ####
#### equivalent to fitting unpenalized domain- and protein-specific intercepts     ####
center_XY <- function(sp){
  muX <- colMeans(sp$X_train)
  muY <- colMeans(sp$Y_train)
  sp$X_train <- sweep(as.matrix(sp$X_train), 2, muX)
  sp$X_val   <- sweep(as.matrix(sp$X_val),   2, muX)
  sp$X_test  <- sweep(as.matrix(sp$X_test),  2, muX)
  sp$Y_train <- sweep(as.matrix(sp$Y_train), 2, muY)
  sp$Y_val   <- sweep(as.matrix(sp$Y_val),   2, muY)
  sp$Y_test  <- sweep(as.matrix(sp$Y_test),  2, muY)
  sp$muX     <- muX
  sp$muY     <- muY
  sp
}

p = ncol(X0) ; K = ncol(Y0)

load(paste('../MTL_transfer_cluster/summary_logtransform_target',tgtidx,'.Rdata',sep=""))

prediction_mse_vec = rep(NA,nrep)
W_opt_debiased_list = list()

for(j in 1:nrep){
  
  set.seed(j)
  
  cat('replication = ',j,'\n')
  
  X0_split_data = center_XY(split_data(X0, Y0, train_prop, val_prop, test_prop))
  
  W_opt_fused = W_opt_fused_list[[j]]
  
  #### Then, debiasing ####
  #### choose lambda_debias by cross validation ####
  
  cv_result = expand.grid(lambda_values) # grid for lam_debias
  cv_result[,2] = NA
  colnames(cv_result) = c('lam_debias','cvrmse')
  
  X_target = X0_split_data$X_train ; Y_target = X0_split_data$Y_train
  
  for(i in 1:length(lambda_values)){ # start cv
    
    lam_debias = cv_result[i,'lam_debias']
    admmm_mtl = MTL_admm(X = X_target, Y = Y_target - X_target %*% W_opt_fused, lambda = lam_debias, rho = 1, max_iter = max_iter_vanilla, tol = tol_vanilla, verbose = T)
    W_opt_bias = admmm_mtl$sol
    W_opt_debiased = W_opt_fused + W_opt_bias
    
    val_rmse = rmse(X0_split_data$X_val %*% W_opt_debiased, X0_split_data$Y_val)
    cv_result[i,'cvrmse'] = val_rmse
    
  } # end debiasing lambda CV
  
  print(cv_result)
  
  ## best lambda_debiased
  lambda_debias = cv_result[which.min(cv_result$cvrmse), 'lam_debias']
  cat('iteration = ',j, ',lambda_debias = ', lambda_debias, '\n')
  
  #### predict with optimal values ####
  X_target = rbind(X0_split_data$X_train, X0_split_data$X_val) ; Y_target = rbind(X0_split_data$Y_train, X0_split_data$Y_val)
  
  admmm_mtl = MTL_admm(X = X_target, Y = Y_target - X_target %*% W_opt_fused, lambda = lambda_debias, rho = 1, max_iter = max_iter_vanilla, tol = tol_vanilla, verbose = T)
  W_opt_bias = admmm_mtl$sol
  W_opt_debiased = W_opt_fused + W_opt_bias
  
  #parameter_list = W_opt_debiased
  W_opt_debiased_list[[j]] = W_opt_debiased
  
  test_predictions = sweep(as.matrix(X0_split_data$X_test %*% W_opt_debiased), 2, X0_split_data$muY, "+")
  Y_test_true = sweep(as.matrix(X0_split_data$Y_test), 2, X0_split_data$muY, "+")
  
  prediction_mse_vec[j] = sse(test_predictions, Y_test_true)
  
  print(prediction_mse_vec)
  cat(100*mean(rowMeans(W_opt_debiased) == 0),'% of the features are zero','\n')
  
  save(W_opt_debiased_list, prediction_mse_vec, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
}
save(W_opt_debiased_list, prediction_mse_vec, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))