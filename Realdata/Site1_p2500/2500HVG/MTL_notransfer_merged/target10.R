library(glmnet)
library(Matrix)
library(dplyr)
source('functions_joint.R')
source('simulationsetting.R')

tgtidx = 10

#### load each domain separately (kept per domain for domain-wise centering) ####
Xlist = vector('list', 12) ; Ylist = vector('list', 12)
for(i in 1:12){
  load(paste('../clr_data_donor_nofiltered_2500hvg',i,'.Rdata',sep=""))
  Xlist[[i]] = as.matrix(X) ; Ylist[[i]] = as.matrix(Y)
}
rm(X,Y,i)

print(dim(Xlist[[tgtidx]])) ; print(dim(Ylist[[tgtidx]]))

p = ncol(Xlist[[1]]) ; K = ncol(Ylist[[1]])

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

parameter_list = list()
prediction_mse_vec = rep(NA,nrep)

for(j in 1:nrep){
  
  cat('replication = ',j,'\n')
  
  #### split each domain separately with the same seed, then pool          ####
  #### the target block is then identical to the MTL(Target) split         ####
  sp = vector('list', 12)
  for(i in 1:12){
    set.seed(j)
    sp[[i]] = center_XY(split_data(Xlist[[i]], Ylist[[i]], train_prop, val_prop, test_prop))
  }
  
  X_tr = do.call(rbind, lapply(sp, function(s) s$X_train))
  Y_tr = do.call(rbind, lapply(sp, function(s) s$Y_train))
  X_va = do.call(rbind, lapply(sp, function(s) s$X_val))
  Y_va = do.call(rbind, lapply(sp, function(s) s$Y_val))
  
  X0__split_data = sp[[tgtidx]]
  
  cv_result = expand.grid(lambda_values)
  cv_result[,2] = NA
  colnames(cv_result) = c('lam','cvrmse')
  
  for(i in 1:nrow(cv_result)){
    cat('lambda = ', lambda_values[i], '\n')
    lam = cv_result[i, 'lam']
    W_opt = MTL_admm(X = X_tr, Y = Y_tr, lambda = lam, rho = 1, max_iter = max_iter_vanilla_merged, tol = tol_vanilla_merged)$sol
    cv_result[i, 'cvrmse'] = rmse(X0__split_data$X_val %*% W_opt, X0__split_data$Y_val)
  }
  
  print(cv_result)
  
  lambda = cv_result[which.min(cv_result$cvrmse), 'lam'] ; cat('iteration = ',j, 'lambda = ', lambda, '\n')
  
  W_opt = MTL_admm(X = rbind(X_tr, X_va), Y = rbind(Y_tr, Y_va), lambda = lambda, rho = 1, max_iter = max_iter_vanilla_merged, tol = tol_vanilla_merged)$sol
  
  parameter_list[[j]] = W_opt
  test_predictions = sweep(as.matrix(X0__split_data$X_test %*% W_opt), 2, X0__split_data$muY, "+")
  Y_test_true = sweep(as.matrix(X0__split_data$Y_test), 2, X0__split_data$muY, "+")
  prediction_mse_vec[j] = sse(test_predictions, Y_test_true)
  
  print(prediction_mse_vec)
  cat(100*mean(rowMeans(W_opt) == 0),'% of the features are zero','\n')
  save( prediction_mse_vec, parameter_list, lambda, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
}
save( prediction_mse_vec, parameter_list, lambda, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))