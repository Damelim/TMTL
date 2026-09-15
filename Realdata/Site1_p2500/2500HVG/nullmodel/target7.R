library(Matrix)
library(dplyr)

#source('../../../../../../functions_joint.R')
source('simulationsetting.R')

tgtidx = 7
load(paste('../clr_data_donor_nofiltered_2500hvg',tgtidx,'.Rdata',sep=""))
X0 = X ; Y0 = Y ; Y0 = as.matrix(Y0) 

rm(X, Y) 

#### donor-wise, training-only centering of the CLR responses ####
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
  
  #### null model : W = 0, i.e. predict the training mean of each protein ####
  W_opt = matrix(0, p, K)
  
  parameter_list[[j]] = W_opt
  test_predictions = sweep(as.matrix(X0_split_data$X_test %*% W_opt), 2, X0_split_data$muY, "+")
  Y_test_true = sweep(as.matrix(X0_split_data$Y_test), 2, X0_split_data$muY, "+")
  prediction_mse_vec[j] = sse(test_predictions, Y_test_true)
  
  print(prediction_mse_vec)
  
  save(parameter_list, prediction_mse_vec, test_sample_size, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
}
save(parameter_list, prediction_mse_vec, test_sample_size, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))