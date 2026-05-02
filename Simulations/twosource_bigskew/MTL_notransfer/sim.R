library(glmnet)
library(Matrix)
library(dplyr)

#### Model 2) MTL with no transfer ####

parameter_list = list()
estimation_msevec = rep(NA, nreplications)
estimation_mse_bytask_list = list()
prediction_msevec = rep(NA, nreplications)
prediction_mse_bytask_list = list()

for(j in 1:nreplications){ # replication loop
  
  set.seed(j)
  cat('replication =', j, '\n')
  
  ## Generate (X,Y)_src, (X,Y)_tgt
  X_target = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) # ; X_target = scale(X_target, center = T, scale = F) #; X_target_center = attr(X_target, "scaled:center")
  X_target_test = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) # iid copy of X_target to evaluate out of sample RMSE.
  Y_target = X_target %*% W_true_target + matrix(rnorm(n_target * K, mean = 0, sd = sigma), n_target, K) #; Y_target = scale(Y_target, center = T, scale = F) ; Y_target_center = attr(Y_target, "scaled:center")
  
  X_source1 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) # ; X_source1 = scale(X_source1, center = T, scale = F) ; X_source1_center = attr(X_source1, "scaled:center")
  Y_source1 = X_source1 %*% W_true_source1 + matrix(rnorm(n_source * K, sd = sigma), n_source, K) # ; Y_source1 = scale(Y_source1, center = T, scale = F) ; Y_source1_center = attr(Y_source1, "scaled:center")
  X_source2 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) # ; X_source2 = scale(X_source2, center = T, scale = F) ; X_source2_center = attr(X_source2, "scaled:center")
  Y_source2 = X_source2 %*% W_true_source2 + matrix(rnorm(n_source * K, sd = sigma), n_source, K) # ; Y_source2 = scale(Y_source2, center = T, scale = F) ; Y_source2_center = attr(Y_source2, "scaled:center")
  
  cv_result = expand.grid(lambda_values)
  cv_result[,2] = NA
  colnames(cv_result) = c('lam','cvrmse')
  cv_rmse_list = rep(NA, cv_folds)
  
  for(i in 1:nrow(cv_result)){ # start cv
    
    lam = cv_result[i,'lam'] 
    
    fold_idx_target = sample(rep(1:cv_folds, length.out = nrow(X_target)))
    
    for(fold in 1:cv_folds){
      train_target_idx = which(fold_idx_target != fold) ; val_target_idx = which(fold_idx_target == fold)
      X_target_tr = X_target[train_target_idx,] ; X_target_val = X_target[val_target_idx,]
      Y_target_tr = Y_target[train_target_idx,] ; Y_target_val = Y_target[val_target_idx,]
      
      W_opt = MTL_admm(X = X_target_tr, Y = Y_target_tr, lambda = lam,
                       rho = 1, max_iter = 1000, tol = 1e-5)$sol
      
      val_rmse = rmse(X_target_val %*% W_opt, Y_target_val)
      cv_rmse_list[fold] = val_rmse
      
    }# end fold loop
      
    avg_cv_rmse = mean(cv_rmse_list)
    cv_result[i,'cvrmse'] = avg_cv_rmse
    
  }# end cv
  
  # Best lambda and fit model with best chosen lambda
  lambda = cv_result[which.min(cv_result$cvrmse), 'lam']
  W_opt = MTL_admm(X = X_target, Y = Y_target, lambda = lambda,
                   rho = 1, max_iter = 1000, tol = 1e-5)$sol
  parameter_list[[j]] = W_opt
  estimation_mse_bytask_list[[j]] = apply(W_true_target - W_opt, 2, function(x) sqrt(mean((x)^2)))
  estimation_msevec[j] = rmse(W_true_target, W_opt)
  prediction_msevec[j] = rmse(X_target_test %*% W_opt, X_target_test %*% W_true_target)
  prediction_mse_bytask_list[[j]] = colMeans((X_target_test %*% W_opt - X_target_test %*% W_true_target)^2) %>% sqrt()
  
  save(parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, 
       file = "summary.Rdata")
  
}# end replication loop

save(parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, 
     file = "summary.Rdata")