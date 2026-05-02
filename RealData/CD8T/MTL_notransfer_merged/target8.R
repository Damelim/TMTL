library(glmnet)
library(Matrix)
library(dplyr)

set.seed(1)

source('~/Library/CloudStorage/GoogleDrive-96limtotoro@gmail.com/My Drive/research_multitask/functions_highdim.R')

load('../clr_data_donor_nofiltered_1000hvg8.Rdata')
X0_ = X ; Y0_ = Y ;  X0_ = as(X0_, "dgCMatrix") ; Y0_ = as(Y0_, "dgCMatrix")
print(dim(X0_)) ; print(dim(Y0_))

#### merge data #####
load('../clr_data_donor_nofiltered_1000hvg1.Rdata')
X0 = X ; Y0 = Y
for(i in 2:8){
  load(paste('../clr_data_donor_nofiltered_1000hvg',i,'.Rdata',sep=""))
  X0 = rbind(X0,X) ; Y0 = rbind(Y0,Y)
}
X0 = as(X0, "dgCMatrix") ; Y0 = as(Y0, "dgCMatrix")

rm(X,Y,i)

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



prediction_mse_vec = rep(NA,nrep)

for(j in 1:nrep){
  
  cat('replication = ',j,'\n')
  
  X0__split_data = split_data(X0_, Y0_, train_prop, val_prop, test_prop, seed = j-1)
  X0_split_data = split_data(X0, Y0, train_prop, val_prop, test_prop, seed = j-1+5)
  
  cv_result = expand.grid(lambda_values)
  cv_result[,2] = NA
  colnames(cv_result) = c('lam','cvrmse')
  
  
  for(i in 1:nrow(cv_result)){
    cat('lambda = ', lambda_values[i], '\n')
    lam = cv_result[i, 'lam']
    W_opt = MTL_admm(X = X0_split_data$X_train, Y = X0_split_data$Y_train, lambda = lam, rho = 1, max_iter = 50, tol = sqrt(K)*1e-3, sparse = T, verbose=T)$sol
    cv_result[i, 'cvrmse'] = rmse(X0__split_data$X_val %*% W_opt, X0__split_data$Y_val)
  }
  
  lambda = cv_result[which.min(cv_result$cvrmse), 'lam'] ; cat('iteration = ',j, 'lambda = ', lambda, '\n')
  
  
  
  W_opt = MTL_admm(X = rbind(X0_split_data$X_train, X0_split_data$X_val), Y = rbind(X0_split_data$Y_train, X0_split_data$Y_val), lambda = lambda, rho = 1, max_iter = 50, tol = sqrt(K)*1e-3, sparse = T, verbose=T)$sol
  parameter_list = W_opt
  test_predictions = X0__split_data$X_test %*% W_opt ; Y_test_true = X0__split_data$Y_test
  prediction_mse_vec[j] = sse(test_predictions, Y_test_true)
  
  print(prediction_mse_vec)
  cat(100*mean(rowMeans(W_opt) == 0),'% of the features are zero','\n')
  save( prediction_mse_vec, parameter_list, lambda, file = "summary_logtransform_target8.Rdata")
}

save( prediction_mse_vec, parameter_list, lambda, file = "summary_logtransform_target8.Rdata")