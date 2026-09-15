library(glmnet)
library(Matrix)
library(dplyr)
#### Model 3) He et al. 2024 (B^{debiased}) ####
#### fused estimator는 ../STL_transfer_cluster_sparse 에서 이미 계산됨. 여기서는 debiasing만 수행.

fusedenv <- new.env()
load('../STL_transfer_cluster_sparse/summary.Rdata', envir = fusedenv)
W_fused_list <- get('parameter_list', envir = fusedenv)   # [[j]] = p x K fused coefficients
# B0_list <- get('B0_list', envir = fusedenv)
# B1_list <- get('B1_list', envir = fusedenv)
# B2_list <- get('B2_list', envir = fusedenv)
# B3_list <- get('B3_list', envir = fusedenv)
lambda1_list <- get('lambda1_list', envir = fusedenv)      # fused 단계에서 선택된 lambda0

nrep_avail <- sum(!sapply(W_fused_list, is.null))
nrep_use   <- min(nreplications, nrep_avail)
cat('fused replications available =', nrep_avail, '/ using', nrep_use, '\n')

parameter_list = list()
estimation_msevec = rep(NA, nreplications)
estimation_mse_bytask_list = list()
prediction_msevec = rep(NA, nreplications)
prediction_mse_bytask_list = list()
lambda_debias_list = list()

for(j in 1:nrep_use){ # replication loop
  
  set.seed(j)
  cat('replication = ', j, '\n')
  
  ## Generate (X,Y)_src, (X,Y)_tgt  -- fused 스크립트와 동일한 순서/호출 (RNG 일치 목적)
  X_target = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p)
  X_target_test = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) # iid copy of X_target to evaluate out of sample RMSE.
  Y_target = X_target %*% W_true_target + matrix(rnorm(n_target * K, mean = 0, sd = sigma), n_target, K)
  
  X_source1 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) 
  Y_source1 = X_source1 %*% W_true_source1 + matrix(rnorm(n_source * K, sd = sigma), n_source, K)  
  X_source2 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p)  
  Y_source2 = X_source2 %*% W_true_source2 + matrix(rnorm(n_source * K, sd = sigma), n_source, K)  
  X_source3 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p)  
  Y_source3 = X_source3 %*% W_true_source3 + matrix(rnorm(n_source * K, sd = sigma), n_source, K)  
  
  W_fused = W_fused_list[[j]]   # p x K, 이미 계산된 fused estimator
  
  singletask_target_rmse_bytask_transfer = rep(NA, K)
  best_lambda_debias_single = rep(NA, K)
  singletask_obtained_coefficients_transfer = matrix(NA, p, K)
  
  for(k in 1:K){ # task loop
    
    if(k %% 5 == 0){
      cat('task', k, '\n')
    }
    
    w_opt_fused = W_fused[,k]
    
    # debiasing
    resid_k = Y_target[,k] - X_target %*% w_opt_fused
    cv_fit <- cv.glmnet(X_target, resid_k, lambda = lambda_values, alpha = 1, nfolds = 5, family = "gaussian")
    model = glmnet(X_target, resid_k, alpha = 1,
                   lambda = cv_fit$lambda.min, family = "gaussian") 
    w_opt_debiased = w_opt_fused + as.numeric(model$beta)
    
    best_lambda_debias_single[k] = cv_fit$lambda.min
    singletask_obtained_coefficients_transfer[,k] = w_opt_debiased
    singletask_target_rmse_bytask_transfer[k] = rmse(X_target_test %*% w_opt_debiased, X_target_test %*% W_true_target[,k]) 
  }# end task loop
  
  parameter_list[[j]] = singletask_obtained_coefficients_transfer
  estimation_msevec[j] = rmse(W_true_target, singletask_obtained_coefficients_transfer)
  estimation_mse_bytask_list[[j]] = apply(W_true_target - singletask_obtained_coefficients_transfer, 2, function(x) sqrt(mean((x)^2)))
  prediction_msevec[j] = sqrt(mean(singletask_target_rmse_bytask_transfer^2))
  prediction_mse_bytask_list[[j]] = singletask_target_rmse_bytask_transfer
  lambda_debias_list[[j]] = best_lambda_debias_single
  
  save(parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda_debias_list,
       file = "summary.Rdata")
  
}# end iteration loop
save(parameter_list, estimation_msevec, estimation_mse_bytask_list, prediction_msevec, prediction_mse_bytask_list, lambda1_list, lambda_debias_list,
     file = "summary.Rdata")