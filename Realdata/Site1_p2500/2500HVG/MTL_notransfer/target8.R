library(Matrix)
library(dplyr)
source('functions_joint.R')
source('simulationsetting.R')

tgtidx = 8

load(paste('../clr_data_donor_nofiltered_2500hvg',tgtidx,'.Rdata',sep=""))
X0 = X ; Y0 = Y ; Y0 = as.matrix(Y0) ; names(which(is.na(colSums(X0))))
hist(as.numeric(X0))
hist(as.numeric(Y0))
print(dim(X0)) ; print(dim(Y0))
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

prediction_mse_vec = rep(NA,nrep)
parameter_list = list()

for(j in 1:nrep){
  
  set.seed(j)
  
  cat('replication = ',j,'\n')
  X0_split_data = center_XY(split_data(X0, Y0, train_prop, val_prop, test_prop))
  
  test_sample_size = nrow(X0_split_data$Y_test)
  
  cv_result = expand.grid(lambda_values)
  cv_result[,2] = NA
  colnames(cv_result) = c('lam','cvrmse')
  
  for(i in 1:nrow(cv_result)){
    cat('lambda = ', lambda_values[i], '\n')
    lam = cv_result[i, 'lam']
    W_opt = MTL_admm(X = X0_split_data$X_train, Y = X0_split_data$Y_train, lambda = lam, rho = 1, max_iter = max_iter_vanilla, tol = tol_vanilla, 
                     verbose=T)$sol  
    cv_result[i, 'cvrmse'] = rmse(X0_split_data$X_val %*% W_opt, X0_split_data$Y_val) # no need to expm1 retransform here because I just need to compare which is the biggest. expm1(.) is monotone.
  }
  
  print(cv_result)
  
  lambda = cv_result[which.min(cv_result$cvrmse), 'lam'] ; cat('iteration = ',j, 'lambda = ', lambda, '\n')
  
  W_opt = MTL_admm(X = rbind(X0_split_data$X_train, X0_split_data$X_val), Y = rbind(X0_split_data$Y_train, X0_split_data$Y_val), lambda = lambda, rho = 1, max_iter = max_iter_vanilla, tol = tol_vanilla, 
                   verbose = T)$sol  
  
  parameter_list[[j]] = W_opt
  test_predictions = sweep(as.matrix(X0_split_data$X_test %*% W_opt), 2, X0_split_data$muY, "+")
  Y_test_true = sweep(as.matrix(X0_split_data$Y_test), 2, X0_split_data$muY, "+")
  prediction_mse_vec[j] = sse(test_predictions, Y_test_true)
  
  print(prediction_mse_vec)
  cat(100*mean(rowMeans(W_opt) == 0),'% of the features are zero','\n')
  
  save(parameter_list, lambda, prediction_mse_vec, test_sample_size, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
}
save(parameter_list, lambda, prediction_mse_vec, test_sample_size, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))