library(glmnet)
library(Matrix)
library(dplyr)

set.seed(1)
rmse = function(x,y){
  sqrt(mean((x-y)^2))
}
sse = function(x,y){
  sum((x-y)^2)
}
tgtidx = 2
load(paste('../clr_data_donor_nofiltered_1000hvg',tgtidx,'.Rdata',sep=""))
X0 = X ; Y0 = Y  

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

load(paste('../STL_transfer_cluster/summary_logtransform_target',tgtidx,'.Rdata',sep=""))
prediction_mse_vec = rep(NA,nrep)

W_opt_debiased_list = list()
#### Model 3) STL B^debiased ####
for(j in 1:nrep){
  
  W_opt_debiased = matrix(NA,p,K)
  cat('replication = ',j,'\n')
  X0_split_data = split_data(X0, Y0, train_prop, val_prop, test_prop, seed = j-1)
  X_target = rbind(X0_split_data$X_train, X0_split_data$X_val) ; Y_target = rbind(X0_split_data$Y_train, X0_split_data$Y_val)
  W_opt_fused = W_opt_fused_list[[j]]
  
  for(k in 1:K){ # task loop
    
    if(k%%10==0){cat('k = ',k,'\n')}
    
    w_opt_fused = W_opt_fused[,k]
    
    # debiasing
    cv_fit <- cv.glmnet(X_target, Y_target[,k] - X_target %*% w_opt_fused, lambda = lambda_values, alpha = 1, nfolds = 3, family = "gaussian")
    model = glmnet(X_target, Y_target[,k] - X_target %*% w_opt_fused, alpha = 1,
                   lambda = cv_fit$lambda.min, family = "gaussian") 
    w_opt_debiased = w_opt_fused + as.numeric(model$beta)
    if(k %% 10 == 0){
      cat('k=',k,'lambda=',cv_fit$lambda.min,'\n')
      print(round(as.numeric(model$beta),3))
    }
    
    W_opt_debiased[,k] = w_opt_debiased
  } # end task loop
  
  W_opt_debiased_list[[j]] = W_opt_debiased
  test_predictions = X0_split_data$X_test %*% W_opt_debiased ; Y_test_true = X0_split_data$Y_test
  prediction_mse_vec[j] = sse(test_predictions, X0_split_data$Y_test)
  print(prediction_mse_vec)
  save(prediction_mse_vec, W_opt_debiased_list, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
}# end iteration loop
save(prediction_mse_vec, W_opt_debiased_list, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))