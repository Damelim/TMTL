library(glmnet)
library(Matrix)
library(dplyr)

#### Model 1) MTL with transfer (debiased) ####

B0_list = list()
B1_list = list()
B2_list = list()
parameter_list = list()
estimation_msevec = rep(NA, nreplications)
estimation_mse_bytask_list = list()
prediction_msevec = rep(NA, nreplications)
prediction_mse_bytask_list = list()

lambda1_list = rep(NA, nreplications) ; lambda2_list = rep(NA, nreplications) ; lambda3_list = rep(NA, nreplications) ; lambda_debias_list = rep(NA, nreplications)

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
  
  X0tX0 = crossprod(X_target,X_target) ; X0tY0 = crossprod(X_target,Y_target)
  X1tX1 = crossprod(X_source1,X_source1) ; X1tY1 = crossprod(X_source1,Y_source1)
  X2tX2 = crossprod(X_source2,X_source2) ; X2tY2 = crossprod(X_source2,Y_source2)
  
  cv_result = expand.grid(lambda_values, lambda_values, lambda_values)
  cv_result[,4] = NA
  colnames(cv_result) = c('lam1','lam2','lam3','cvrmse')
  cv_rmse_list = rep(NA, cv_folds)

  for(i in 1:nrow(cv_result)){ # start cv
    
    if(i %% 50 == 0){ cat('i = ', i, '\n') }
    
    lam1 = cv_result[i,'lam1'] ; lam2 = cv_result[i,'lam2'] ; lam3 = cv_result[i,'lam3']
    
    fold_idx_source1 = sample(rep(1:cv_folds, length.out = nrow(X_source1)))
    fold_idx_source2 = sample(rep(1:cv_folds, length.out = nrow(X_source2)))
    fold_idx_target = sample(rep(1:cv_folds, length.out = nrow(X_target)))
    
    for(fold in 1:cv_folds){
      
      train_source1_idx = which(fold_idx_source1 != fold) ; val_source1_idx = which(fold_idx_source1 == fold)
      train_source2_idx = which(fold_idx_source2 != fold) ; val_source2_idx = which(fold_idx_source2 == fold)
      train_target_idx = which(fold_idx_target != fold) ; val_target_idx = which(fold_idx_target == fold)
      
      X_source1_tr = X_source1[train_source1_idx,] ; X_source1_val = X_source1[val_source1_idx,]
      Y_source1_tr = Y_source1[train_source1_idx,] ; Y_source1_val = Y_source1[val_source1_idx,]
      X_source2_tr = X_source2[train_source2_idx,] ; X_source2_val = X_source2[val_source2_idx,]
      Y_source2_tr = Y_source2[train_source2_idx,] ; Y_source2_val = Y_source2[val_source2_idx,]
      X_target_tr = X_target[train_target_idx,] ; X_target_val = X_target[val_target_idx,]
      Y_target_tr = Y_target[train_target_idx,] ; Y_target_val = Y_target[val_target_idx,]
      
      X0tX0_tr = crossprod(X_target_tr, X_target_tr) ; X0tY0_tr = crossprod(X_target_tr, Y_target_tr)
      X1tX1_tr = crossprod(X_source1_tr, X_source1_tr) ; X1tY1_tr = crossprod(X_source1_tr, Y_source1_tr)
      X2tX2_tr = crossprod(X_source2_tr, X_source2_tr) ; X2tY2_tr = crossprod(X_source2_tr, Y_source2_tr)
      
      if(i %% 50 == 0){
        admmm = admm_two_source(X0 = X_target_tr, Y0 = Y_target_tr, X1 = X_source1_tr, Y1 = Y_source1_tr, X2 = X_source2_tr, Y2 = Y_source2_tr,
                                X0tX0_tr, X0tY0_tr, X1tX1_tr, X1tY1_tr, X2tX2_tr, X2tY2_tr,
                                lambda0 = lam1, lambda1 = lam2, lambda2 = lam3, 
                                rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1,
                                mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = T)
      }else{
        admmm = admm_two_source(X0 = X_target_tr, Y0 = Y_target_tr, X1 = X_source1_tr, Y1 = Y_source1_tr, X2 = X_source2_tr, Y2 = Y_source2_tr,
                                X0tX0_tr, X0tY0_tr, X1tX1_tr, X1tY1_tr, X2tX2_tr, X2tY2_tr,
                                lambda0 = lam1, lambda1 = lam2, lambda2 = lam3, 
                                rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1,
                                mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)
      }
      W_opt_fused = n_source/N*(admmm$Gamma1 + admmm$Gamma2) + n_target/N*admmm$Gamma00    
      
      #validation rmse
      val_rmse = rmse(X_target_val %*% W_opt_fused, Y_target_val)
      cv_rmse_list[fold] = val_rmse
      
    }# end fold loop
    
    avg_cv_rmse = mean(cv_rmse_list)
    cv_result[i,'cvrmse'] = avg_cv_rmse
    
  } # end cv
      
  # Best lambda values
  lambda1 = cv_result[which.min(cv_result$cvrmse), 'lam1'] ; lambda2 = cv_result[which.min(cv_result$cvrmse), 'lam2'] ; lambda3 = cv_result[which.min(cv_result$cvrmse), 'lam3']
  cat('lambda1 = ', lambda1, ', lambda2 = ', lambda2, 'lambda3 = ', lambda3, '\n')
  lambda1_list[j] = lambda1 ; lambda2_list[j] = lambda2 ; lambda3_list[j] = lambda3
  
  #### choose lambda_debias by cross validation ####
  ## first, fit W_opt_fused.
  admmm = admm_two_source(X0 = X_target, Y0 = Y_target, X1 = X_source1, Y1 = Y_source1, X2 = X_source2, Y2 = Y_source2,
                          X0tX0, X0tY0, X1tX1, X1tY1, X2tX2, X2tY2,
                          lambda0 = lambda1, lambda1 = lambda2, lambda2 = lambda3, 
                          rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1,
                          mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                          max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)
  B0_list[[j]] = admmm$Gamma00 ; B1_list[[j]] = admmm$Gamma1 ; B2_list[[j]] = admmm$Gamma2
  W_opt_fused = n_source/N*(admmm$Gamma1 + admmm$Gamma2) + n_target/N*admmm$Gamma00
  
  # Now debiasing
  cv_result = expand.grid(lambda_values) # grid for lam_debias
  cv_result[,2] = NA
  colnames(cv_result) = c('lam_debias','cvrmse')
  cv_rmse_list = rep(NA, cv_folds)
  
  for(i in 1:length(lambda_values)){
    
    lam_debias = cv_result[i,'lam_debias']
   
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
  print(W_opt_bias)
  W_opt_debiased = W_opt_fused + W_opt_bias
  
  parameter_list[[j]] = W_opt_debiased
  estimation_msevec[j] = rmse(W_true_target, W_opt_debiased)
  estimation_mse_bytask_list[[j]] = apply(W_true_target - W_opt_debiased, 2, function(x) sqrt(mean((x)^2)))
    
  prediction_msevec[j] = rmse(X_target_test %*% W_opt_debiased, X_target_test %*% W_true_target)
  prediction_mse_bytask_list[[j]] = colMeans((X_target_test %*% W_opt_debiased - X_target_test %*% W_true_target)^2) %>% sqrt()
  
  save(B0_list, B1_list, B2_list, parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda2_list, lambda3_list, lambda_debias_list,
       file = "summary.Rdata")
}

save(B0_list, B1_list, B2_list, parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda2_list, lambda3_list, lambda_debias_list,
     file = "summary.Rdata")