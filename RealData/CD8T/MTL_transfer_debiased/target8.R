library(glmnet)
library(Matrix)
library(dplyr)

set.seed(1)

source('~/Library/CloudStorage/GoogleDrive-96limtotoro@gmail.com/My Drive/research_multitask/functions.R')

tgtidx = 8
load(paste('../clr_data_donor_nofiltered_1000hvg',tgtidx,'.Rdata',sep=""))
X0 = X ; Y0 = Y ; Y0 = as.matrix(Y0) #Y0 = as.matrix(log1p(Y0)) #; X0 = as(X0, "dgCMatrix") ; Y0 = as(Y0, "dgCMatrix")

rm(X, Y) 

split_data <- function(X, Y, train_prop = 0.6, val_prop = 0.2, test_prop = 0.2, seed) {
  set.seed(seed)
  n <- nrow(X)
  idx_all <- sample(n)   # Shuffle indices
  
  n_train <- floor(train_prop * n)
  n_val <- floor(val_prop * n)
  n_test <- n - n_train - n_val
  
  idx_train <- idx_all[1:n_train]
  idx_val <- idx_all[(n_train + 1):(n_train + n_val)]
  idx_test <- idx_all[(n_train + n_val + 1):n]
  
  return(list(
    X_train = X[idx_train, , drop = FALSE],
    Y_train = Y[idx_train, , drop = FALSE],
    X_val   = X[idx_val, , drop = FALSE],
    Y_val   = Y[idx_val, , drop = FALSE],
    X_test  = X[idx_test, , drop = FALSE],
    Y_test  = Y[idx_test, , drop = FALSE]
  ))
}

p = ncol(X0) ; K = ncol(Y0)

load(paste('../MTL_transfer/summary_logtransform_target',tgtidx,'.Rdata',sep=""))

prediction_mse_vec = rep(NA,nrep)
W_opt_debiased_list = list()

for(j in 1:nrep){
  
  cat('replication = ',j,'\n')
  
  X0_split_data = split_data(X0, Y0, train_prop, val_prop, test_prop, seed = j-1)
  
  W_opt_fused = W_opt_fused_list[[j]]
  
  #### Then, debiasing ####
  #### choose lambda_debias by cross validation ####
  
  cv_result = expand.grid(lambda_values) # grid for lam_debias
  cv_result[,2] = NA
  colnames(cv_result) = c('lam_debias','cvrmse')
  
  X_target = X0_split_data$X_train ; Y_target = X0_split_data$Y_train
  
  for(i in 1:length(lambda_values)){ # start cv
    
    lam_debias = cv_result[i,'lam_debias']
    admmm_mtl = MTL_admm(X = X_target, Y = Y_target - X_target %*% W_opt_fused, lambda = lam_debias, rho = 1, max_iter = 100, tol = sqrt(K)*5e-5, verbose = T)
    W_opt_bias = admmm_mtl$sol
    W_opt_debiased = W_opt_fused + W_opt_bias
    
    val_rmse = rmse(X0_split_data$X_val %*% W_opt_debiased, X0_split_data$Y_val)
    cv_result[i,'cvrmse'] = val_rmse
    
  } # end debiasing lambda CV
  
  ## best lambda_debiased
  lambda_debias = cv_result[which.min(cv_result$cvrmse), 'lam_debias']
  cat('iteration = ',j, ',lambda_debias = ', lambda_debias, '\n')
  
  #### predict with optimal values ####
  X_target = rbind(X0_split_data$X_train, X0_split_data$X_val) ; Y_target = rbind(X0_split_data$Y_train, X0_split_data$Y_val)
  
  admmm_mtl = MTL_admm(X = X_target, Y = Y_target - X_target %*% W_opt_fused, lambda = lambda_debias, rho = 1, max_iter = 100, tol = sqrt(K)*5e-5, verbose = T)
  W_opt_bias = admmm_mtl$sol
  W_opt_debiased = W_opt_fused + W_opt_bias
  
  #parameter_list = W_opt_debiased
  W_opt_debiased_list[[j]] = W_opt_debiased
  
  test_predictions = X0_split_data$X_test %*% W_opt_debiased ; Y_test_true = X0_split_data$Y_test
  
  prediction_mse_vec[j] = sse(test_predictions, X0_split_data$Y_test)
  
  print(prediction_mse_vec)
  
  save(W_opt_debiased_list, prediction_mse_vec, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
}
save(W_opt_debiased_list, prediction_mse_vec, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))