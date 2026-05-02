library(glmnet)
library(Matrix)
library(dplyr)

#### Model 1) MTL with transfer (debiasing) ####

B0_list = list()
B1_list = list()
parameter_list = list()
estimation_msevec = rep(NA, nreplications)
estimation_mse_bytask_list = list()
prediction_msevec = rep(NA, nreplications)
prediction_mse_bytask_list = list()

lambda1_list = rep(NA, nreplications) ; lambda2_list = rep(NA, nreplications) ; lambda_debias_list = rep(NA, nreplications)

for(j in 1:nreplications){ # replication loop
  
  set.seed(j)
  cat('replication =', j, '\n')
  
  ## Generate (X,Y)_src, (X,Y)_tgt
  X_target = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) #; X_target = scale(X_target, center = T, scale = F) #; X_target_center = attr(X_target, "scaled:center")
  X_target_test = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) # iid copy of X_target to evaluate out of sample RMSE.
  Y_target = X_target %*% W_true_target + matrix(rnorm(n_target * K, mean = 0, sd = sigma), n_target, K) #; Y_target = scale(Y_target, center = T, scale = F) ; Y_target_center = attr(Y_target, "scaled:center")

  X_source = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) #; X_source = scale(X_source, center = T, scale = F) ; X_source_center = attr(X_source, "scaled:center")
  Y_source = X_source %*% W_true_source + matrix(rnorm(n_source * K, sd = sigma), n_source, K) #; Y_source = scale(Y_source, center = T, scale = F) ; Y_source_center = attr(Y_source, "scaled:center")
  
  X0tX0 = crossprod(X_target,X_target) ; X0tY0 = crossprod(X_target,Y_target)
  X1tX1 = crossprod(X_source,X_source) ; X1tY1 = crossprod(X_source,Y_source)
  
  
  cv_result = expand.grid(lambda_values, lambda_values)
  cv_result[,3] = NA
  colnames(cv_result) = c('lam1','lam2','cvrmse')
  cv_rmse_list = rep(NA, cv_folds)

  for(i in 1:nrow(cv_result)){ # start cv
    
    if(i %% 25 == 0){
      cat('i = ', i, '\n')
    }
    
    lam1 = cv_result[i,'lam1'] ; lam2 = cv_result[i,'lam2']  
    
    fold_idx_source = sample(rep(1:cv_folds, length.out = nrow(X_source))) # because n_s1 = n_s2. !!!
    fold_idx_target = sample(rep(1:cv_folds, length.out = nrow(X_target)))
    
    for(fold in 1:cv_folds){
      
      train_source_idx = which(fold_idx_source != fold) ; val_source_idx = which(fold_idx_source == fold)
      train_target_idx = which(fold_idx_target != fold) ; val_target_idx = which(fold_idx_target == fold)
      
      X_source_tr = X_source[train_source_idx,] ; X_source_val = X_source[val_source_idx,]
      Y_source_tr = Y_source[train_source_idx,] ; Y_source_val = Y_source[val_source_idx,]
      
      X_target_tr = X_target[train_target_idx,] ; X_target_val = X_target[val_target_idx,]
      Y_target_tr = Y_target[train_target_idx,] ; Y_target_val = Y_target[val_target_idx,]
      
      X0tX0_tr = crossprod(X_target_tr, X_target_tr) ; X0tY0_tr = crossprod(X_target_tr, Y_target_tr)
      X1tX1_tr = crossprod(X_source_tr, X_source_tr) ; X1tY1_tr = crossprod(X_source_tr, Y_source_tr)
      
      if(i %% 25 != 0){
        admmm = admm_multitask(X0 = X_target_tr, Y0 = Y_target_tr, X1 = X_source_tr, Y1 = Y_source_tr, 
                               X0tX0_tr, X0tY0_tr, X1tX1_tr, X1tY1_tr,
                               lambda0 = lam1, lambda1 = lam2, 
                               rho00 = 1, rho01 = 1, rho1 = 1, 
                               mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                               max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)
      }else{
        admmm = admm_multitask(X0 = X_target_tr, Y0 = Y_target_tr, X1 = X_source_tr, Y1 = Y_source_tr, 
                               X0tX0_tr, X0tY0_tr, X1tX1_tr, X1tY1_tr,
                               lambda0 = lam1, lambda1 = lam2, 
                               rho00 = 1, rho01 = 1, rho1 = 1, 
                               mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                               max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = T)
      }
      W_opt_fused = n_source/N*admmm$Gamma1 + n_target/N*admmm$Gamma00
      
      #validation rmse
      val_rmse = rmse(X_target_val %*% W_opt_fused, Y_target_val)
      cv_rmse_list[fold] = val_rmse
    }# end fold loop
    
    avg_cv_rmse = mean(cv_rmse_list)
    cv_result[i,'cvrmse'] = avg_cv_rmse
    
  } # end cv
  
  # Best lambda1, lambda2.
  lambda1 = cv_result[which.min(cv_result$cvrmse), 'lam1'] ; lambda2 = cv_result[which.min(cv_result$cvrmse), 'lam2'] 
  cat('lambda1 = ', lambda1, ', lambda2 = ', lambda2,'\n')
  lambda1_list[j] = lambda1 ; lambda2_list[j] = lambda2
  
  
  ## first, fit W_opt_fused.
  admmm = admm_multitask(X0 = X_target, Y0 = Y_target, X1 = X_source, Y1 = Y_source, 
                         X0tX0, X0tY0, X1tX1, X1tY1,
                         lambda0 = lambda1, lambda1 = lambda2, 
                         rho00 = 1, rho01 = 1, rho1 = 1, 
                         mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                         max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)
  
  W_opt_fused = n_source/N*admmm$Gamma1 + n_target/N*admmm$Gamma00
  
  B0_list[[j]] = admmm$Gamma00 ; B1_list[[j]] = admmm$Gamma1
  
  
  
  
  cv_result = expand.grid(lambda_values) # grid for lam_debias
  cv_result[,2] = NA
  colnames(cv_result) = c('lam_debias','cvrmse')
  cv_rmse_list = rep(NA, cv_folds)
  
  
  for(i in 1:length(lambda_values)){ # start cv
    
    lam_debias = cv_result[i,'lam_debias']
    
    fold_idx_source = sample(rep(1:cv_folds, length.out = nrow(X_source))) # because n_s1 = n_s2. !!!
    fold_idx_target = sample(rep(1:cv_folds, length.out = nrow(X_target)))
    
    for(fold in 1:cv_folds){
      
      train_target_idx = which(fold_idx_target != fold) ; val_target_idx = which(fold_idx_target == fold)
      X_target_tr = X_target[train_target_idx,] ; X_target_val = X_target[val_target_idx,]
      Y_target_tr = Y_target[train_target_idx,] ; Y_target_val = Y_target[val_target_idx,]
      
      X0tX0_tr = crossprod(X_target_tr, X_target_tr) ; X0tY0_tr = crossprod(X_target_tr, Y_target_tr)
      
      admmm_mtl = MTL_admm(X = X_target_tr, Y = Y_target_tr - X_target_tr %*% W_opt_fused, lambda = lam_debias, 
                           rho = 1, max_iter = 1000, 
                           tol = 1e-5)
      W_opt_bias = admmm_mtl$sol
      W_opt_debiased = W_opt_fused + W_opt_bias
      
      #validation rmse
      val_rmse = rmse(X_target_val %*% W_opt_debiased, Y_target_val)
      cv_rmse_list[fold] = val_rmse
    }# end fold loop
    
    avg_cv_rmse = mean(cv_rmse_list)
    cv_result[i,'cvrmse'] = avg_cv_rmse

    
  } # end debiasing lambda CV
  
  
  ## best lambda_debiased
  lambda_debias = cv_result[which.min(cv_result$cvrmse), 'lam_debias']
  cat('lambda_debias = ', lambda_debias, '\n')
  lambda_debias_list[j] = lambda_debias
  
  
  #### predict with optimal values ####
  admmm_mtl = MTL_admm(X = X_target, Y = Y_target - X_target %*% W_opt_fused, lambda = lambda_debias,
                       rho = 1, max_iter = 1000, 
                       tol = 1e-5)
  W_opt_bias = admmm_mtl$sol
  W_opt_debiased = W_opt_fused + W_opt_bias
  print(W_opt_bias)
  
  parameter_list[[j]] = W_opt_debiased
  estimation_msevec[j] = rmse(W_true_target, W_opt_debiased)
  estimation_mse_bytask_list[[j]] = apply(W_true_target - W_opt_debiased, 2, function(x) sqrt(mean((x)^2)))
    
  prediction_msevec[j] = rmse(X_target_test %*% W_opt_debiased, X_target_test %*% W_true_target)
  prediction_mse_bytask_list[[j]] = colMeans((X_target_test %*% W_opt_debiased - X_target_test %*% W_true_target)^2) %>% sqrt()
  
  save(B0_list, B1_list, parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda2_list, lambda_debias_list,
       file = "summary.Rdata")
}
save(B0_list, B1_list, parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda2_list, lambda_debias_list,
     file = "summary.Rdata")