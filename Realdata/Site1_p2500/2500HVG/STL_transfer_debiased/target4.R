library(glmnet)
library(Matrix)
library(dplyr)
#source('../../../../../../functions_joint.R')
source('simulationsetting.R')

tgtidx = 4

load(paste('../clr_data_donor_nofiltered_500hvg',tgtidx,'.Rdata',sep=""))
X0 = X ; Y0 = as.matrix(Y)  
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

load(paste('../STL_transfer_cluster/summary_logtransform_target',tgtidx,'.Rdata',sep=""))

prediction_mse_vec = rep(NA,nrep)
W_opt_debiased_list = list()

#### Model 3) STL B^debiased ####
for(j in 1:nrep){
  
  set.seed(j)
  
  W_opt_debiased = matrix(NA,p,K)
  cat('replication = ',j,'\n')
  
  
  X0_split_data = center_XY(split_data(X0, Y0, train_prop, val_prop, test_prop))
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
  test_predictions = sweep(as.matrix(X0_split_data$X_test %*% W_opt_debiased), 2, X0_split_data$muY, "+")
  Y_test_true = sweep(as.matrix(X0_split_data$Y_test), 2, X0_split_data$muY, "+")
  prediction_mse_vec[j] = sse(test_predictions, Y_test_true)
  print(prediction_mse_vec)
  save(prediction_mse_vec, W_opt_debiased_list, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
}# end iteration loop
save(prediction_mse_vec, W_opt_debiased_list, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))