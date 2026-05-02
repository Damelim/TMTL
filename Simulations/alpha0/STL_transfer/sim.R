library(glmnet)
library(Matrix)
library(dplyr)

#### Model 3) He et al. 2024 (B^{fused}) ####

B0_list = list()
B1_list = list()
parameter_list = list()
estimation_msevec = rep(NA, nreplications)
estimation_mse_bytask_list = list()
prediction_msevec = rep(NA, nreplications)
prediction_mse_bytask_list = list()
lambda1_list = list() ; lambda2_list = list()

for(j in 1:nreplications){ # replication loop
  
  set.seed(j)
  cat('replication = ', j, '\n')

  ## Generate (X,Y)_src, (X,Y)_tgt
  X_target = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) #; X_target = scale(X_target, center = T, scale = F) #; X_target_center = attr(X_target, "scaled:center")
  X_target_test = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) # iid copy of X_target to evaluate out of sample RMSE.
  Y_target = X_target %*% W_true_target + matrix(rnorm(n_target * K, mean = 0, sd = sigma), n_target, K) #; Y_target = scale(Y_target, center = T, scale = F) ; Y_target_center = attr(Y_target, "scaled:center")
  
  # Now, p(X_s) = p(X_t). Deal with COVARIATE SHIFT LATER!
  X_source = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) #; X_source = scale(X_source, center = T, scale = F) ; X_source_center = attr(X_source, "scaled:center")
  Y_source = X_source %*% W_true_source + matrix(rnorm(n_source * K, sd = sigma), n_source, K) #; Y_source = scale(Y_source, center = T, scale = F) ; Y_source_center = attr(Y_source, "scaled:center")
  
  X0tX0 = crossprod(X_target,X_target) ; X1tX1 = crossprod(X_source,X_source)
  
  singletask_target_rmse_bytask_transfer = rep(NA, K)
  best_lambda1_single = rep(NA, K)
  best_lambda2_single = rep(NA, K)
  B0_coef = matrix(NA, p, K) ; B1_coef = matrix(NA, p, K)
  singletask_obtained_coefficients_transfer = matrix(NA, p, K)
  
  for(k in 1:K){ # task loop
    
    if(k %% 5 == 0){
      cat('task', k, '\n')
    }
    X0ty0 = crossprod(X_target,Y_target[,k]) ; X1ty1 = crossprod(X_source,Y_source[,k])

    cv_result = expand.grid(lambda_values, lambda_values)
    cv_result[,3] = NA
    colnames(cv_result) = c('lam1','lam2','cvrmse')
    cv_rmse_list = rep(NA, cv_folds)

    for(i in 1:nrow(cv_result)){ # cv loop 

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
        
        X0tX0_tr = crossprod(X_target_tr, X_target_tr) ; X0ty0_tr = crossprod(X_target_tr, Y_target_tr[,k])
        X1tX1_tr = crossprod(X_source_tr, X_source_tr) ; X1ty1_tr = crossprod(X_source_tr, Y_source_tr[,k])
        
        admmm = admm_stl(X0 = X_target_tr, y0 = Y_target_tr[,k], X1 = X_source_tr, y1 = Y_source_tr[,k], 
                         X0tX0_tr, X0ty0_tr, X1tX1_tr, X1ty1_tr,
                         lambda0 = lam1, lambda1 = lam2,
                         rho00 = 1, rho01 = 1, rho1 = 1, 
                         mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                         max_iter = 500, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)
        w_opt_fused = n_source/N*admmm$Gamma1 + n_target/N*admmm$Gamma00 
        
        val_rmse = rmse(X_target_val %*% w_opt_fused, Y_target_val[,k])
        cv_rmse_list[fold] = val_rmse
      }# end fold loop
      
      avg_cv_rmse = mean(cv_rmse_list)
      cv_result[i,'cvrmse'] = avg_cv_rmse
      
    }# end cv

    # best chosen lambda1 and lambda2
    lambda1 = cv_result[which.min(cv_result$cvrmse), 'lam1'] ; lambda2 = cv_result[which.min(cv_result$cvrmse), 'lam2'] 
    best_lambda1_single[k] = lambda1 ; best_lambda2_single[k] = lambda2
    
    ## first, fused estimator.
    admmm = admm_stl(X0 = X_target, y0 = Y_target[,k], X1 = X_source, y1 = Y_source[,k], 
                     X0tX0, X0ty0, X1tX1, X1ty1,
                     lambda0 = lambda1, lambda1 = lambda2,
                     rho00 = 1, rho01 = 1, rho1 = 1, 
                     mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                     max_iter = 500, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)
    w_opt_fused = n_source/N*admmm$Gamma1 + n_target/N*admmm$Gamma00 
    
    B0_coef[,k] = admmm$Gamma00 ; B1_coef[,k] = admmm$Gamma1
    singletask_obtained_coefficients_transfer[,k] = w_opt_fused
    singletask_target_rmse_bytask_transfer[k] = rmse(X_target_test %*% w_opt_fused, X_target_test %*% W_true_target[,k]) 
  }# end task loop
  
  B0_list[[j]] = B0_coef ; B1_list[[j]] = B1_coef
  parameter_list[[j]] = singletask_obtained_coefficients_transfer
  estimation_msevec[j] = rmse(W_true_target, singletask_obtained_coefficients_transfer)
  estimation_mse_bytask_list[[j]] = apply(W_true_target - singletask_obtained_coefficients_transfer, 2, function(x) sqrt(mean((x)^2)))
  prediction_msevec[j] = sqrt(mean(singletask_target_rmse_bytask_transfer^2))
  prediction_mse_bytask_list[[j]] = singletask_target_rmse_bytask_transfer
  lambda1_list[[j]] = best_lambda1_single ; lambda2_list[[j]] = best_lambda2_single
  
  save(B0_list, B1_list, parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda2_list,
       file = "summary.Rdata")
  
}# end iteration loop

save(B0_list, B1_list, parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda2_list,
     file = "summary.Rdata")