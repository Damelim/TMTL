library(glmnet)
library(Matrix)
library(dplyr)

#### Model 3) He et al. 2024 (B^{debiased}) ####

B0_list = list()
B1_list = list()
B2_list = list()
parameter_list = list()
estimation_msevec = rep(NA, nreplications)
estimation_mse_bytask_list = list()
prediction_msevec = rep(NA, nreplications)
prediction_mse_bytask_list = list()
lambda1_list = list() ; lambda2_list = list() ; lambda3_list = list()

for(j in 1:nreplications){ # replication loop
  
  set.seed(j)
  cat('replication = ', j, '\n')
  
  ## Generate (X,Y)_src, (X,Y)_tgt
  X_target = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) # ; X_target = scale(X_target, center = T, scale = F) #; X_target_center = attr(X_target, "scaled:center")
  X_target_test = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) # iid copy of X_target to evaluate out of sample RMSE.
  Y_target = X_target %*% W_true_target + matrix(rnorm(n_target * K, mean = 0, sd = sigma), n_target, K) #; Y_target = scale(Y_target, center = T, scale = F) ; Y_target_center = attr(Y_target, "scaled:center")
  
  X_source1 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) # ; X_source1 = scale(X_source1, center = T, scale = F) ; X_source1_center = attr(X_source1, "scaled:center")
  Y_source1 = X_source1 %*% W_true_source1 + matrix(rnorm(n_source * K, sd = sigma), n_source, K) # ; Y_source1 = scale(Y_source1, center = T, scale = F) ; Y_source1_center = attr(Y_source1, "scaled:center")
  X_source2 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) # ; X_source2 = scale(X_source2, center = T, scale = F) ; X_source2_center = attr(X_source2, "scaled:center")
  Y_source2 = X_source2 %*% W_true_source2 + matrix(rnorm(n_source * K, sd = sigma), n_source, K) # ; Y_source2 = scale(Y_source2, center = T, scale = F) ; Y_source2_center = attr(Y_source2, "scaled:center")
  
  X0tX0 = crossprod(X_target,X_target) ; X1tX1 = crossprod(X_source1,X_source1) ; X2tX2 = crossprod(X_source2,X_source2)
  
  singletask_target_rmse_bytask_transfer = rep(NA, K)
  best_lambda1_single = rep(NA, K)
  best_lambda2_single = rep(NA, K)
  best_lambda3_single = rep(NA, K)
  B0_coef = matrix(NA, p, K) ; B1_coef = matrix(NA, p, K) ; B2_coef = matrix(NA, p, K)
  singletask_obtained_coefficients_transfer = matrix(NA, p, K)
  
  for(k in 1:K){ # task loop
    
    if(k %% 5 == 0){
      cat('task', k, '\n')
    }
    X0ty0 = crossprod(X_target,Y_target[,k]) ; X1ty1 = crossprod(X_source1,Y_source1[,k]) ; X2ty2 = crossprod(X_source2,Y_source2[,k])
    
    cv_result = expand.grid(lambda_values, lambda_values, lambda_values)
    cv_result[,4] = NA
    colnames(cv_result) = c('lam1','lam2','lam3','cvrmse')
    cv_rmse_list = rep(NA, cv_folds)

    for(i in 1:nrow(cv_result)){ # cv loop 

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
        
        X0tX0_tr = crossprod(X_target_tr, X_target_tr) ; X0ty0_tr = crossprod(X_target_tr, Y_target_tr[,k])
        X1tX1_tr = crossprod(X_source1_tr, X_source1_tr) ; X1ty1_tr = crossprod(X_source1_tr, Y_source1_tr[,k])
        X2tX2_tr = crossprod(X_source2_tr, X_source2_tr) ; X2ty2_tr = crossprod(X_source2_tr, Y_source2_tr[,k])
        
        admmm = admm_stl_twosource(X0 = X_target_tr, y0 = Y_target_tr[,k], X1 = X_source1_tr, y1 = Y_source1_tr[,k], X2 = X_source2_tr, y2 = Y_source2_tr[,k],
                                   X0tX0_tr, X0ty0_tr, X1tX1_tr, X1ty1_tr, X2tX2_tr, X2ty2_tr,
                                   lambda0 = lam1, lambda1 = lam2, lambda2 = lam3,
                                   rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1,
                                   mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                   max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)
        w_opt_fused = n_source/N*(admmm$Gamma1+admmm$Gamma2) + n_target/N*admmm$Gamma00 
        
        val_rmse = rmse(X_target_val %*% w_opt_fused, Y_target_val[,k])
        cv_rmse_list[fold] = val_rmse
      }# end fold loop
      
      avg_cv_rmse = mean(cv_rmse_list)
      cv_result[i,'cvrmse'] = avg_cv_rmse
      
    }# end cv
    
    # best chosen lambda values
    lambda1 = cv_result[which.min(cv_result$cvrmse), 'lam1'] ; lambda2 = cv_result[which.min(cv_result$cvrmse), 'lam2'] ; lambda3 = cv_result[which.min(cv_result$cvrmse), 'lam3'] 
    best_lambda1_single[k] = lambda1 ; best_lambda2_single[k] = lambda2 ; best_lambda3_single[k] = lambda3
    
    admmm = admm_stl_twosource(X0 = X_target, y0 = Y_target[,k], X1 = X_source1, y1 = Y_source1[,k], X2 = X_source2, y2 = Y_source2[,k],
                               X0tX0, X0ty0, X1tX1, X1ty1, X2tX2, X2ty2,
                               lambda0 = lambda1, lambda1 = lambda2, lambda2 = lambda3,
                               rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1,
                               mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                               max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)
    w_opt_fused = n_source/N*(admmm$Gamma1+admmm$Gamma2) + n_target/N*admmm$Gamma00 
    
    
    # debiasing
    cv_fit <- cv.glmnet(X_target, Y_target[,k] - X_target %*% w_opt_fused, lambda = lambda_values, alpha = 1, nfolds = 5, family = "gaussian")
    model = glmnet(X_target, Y_target[,k] - X_target %*% w_opt_fused, alpha = 1,
                   lambda = cv_fit$lambda.min, family = "gaussian") 
    w_opt_debiased = w_opt_fused + as.numeric(model$beta)
    
    B0_coef[,k] = admmm$Gamma00 ; B1_coef[,k] = admmm$Gamma1 ; B2_coef[,k] = admmm$Gamma2
    singletask_obtained_coefficients_transfer[,k] = w_opt_debiased
    singletask_target_rmse_bytask_transfer[k] = rmse(X_target_test %*% w_opt_debiased, X_target_test %*% W_true_target[,k]) 
  }# end task loop
  
  B0_list[[j]] = B0_coef ; B1_list[[j]] = B1_coef ; B2_list[[j]] = B2_coef
  parameter_list[[j]] = singletask_obtained_coefficients_transfer
  estimation_msevec[j] = rmse(W_true_target, singletask_obtained_coefficients_transfer)
  estimation_mse_bytask_list[[j]] = apply(W_true_target - singletask_obtained_coefficients_transfer, 2, function(x) sqrt(mean((x)^2)))
  prediction_msevec[j] = sqrt(mean(singletask_target_rmse_bytask_transfer^2))
  prediction_mse_bytask_list[[j]] = singletask_target_rmse_bytask_transfer
  lambda1_list[[j]] = best_lambda1_single ; lambda2_list[[j]] = best_lambda2_single ; lambda3_list[[j]] = best_lambda3_single
  
  save(B0_list, B1_list, B2_list, parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda2_list, lambda3_list,
       file = "summary.Rdata")
  
}# end iteration loop

save(B0_list, B1_list, B2_list, parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda2_list, lambda3_list,
     file = "summary.Rdata")