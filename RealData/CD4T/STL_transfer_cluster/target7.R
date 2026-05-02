library(glmnet)
library(Matrix)
library(dplyr)

set.seed(1)

source('~/research_multitask/functions_stl_cluster.R')

load("~/research_multitask/dataset/clr_data_donor_nofiltered_1000hvg1.Rdata")
X1 = X ; Y1 = Y %>% as.matrix() ; #Y1 = as.matrix(log1p(Y1)) #; X1 = as(X1, "dgCMatrix") ; Y1 = as(Y1, "dgCMatrix")

load("~/research_multitask/dataset/clr_data_donor_nofiltered_1000hvg2.Rdata")
X2 = X ; Y2 = Y %>% as.matrix(); #Y2 = as.matrix(log1p(Y2)) #; X2 = as(X2, "dgCMatrix") ; Y2 = as(Y2, "dgCMatrix")

load("~/research_multitask/dataset/clr_data_donor_nofiltered_1000hvg3.Rdata")
X3 = X ; Y3 = Y %>% as.matrix() ; #Y3 = as.matrix(log1p(Y3)) #; X3 = as(X3, "dgCMatrix") ; Y3 = as(Y3, "dgCMatrix")

load("~/research_multitask/dataset/clr_data_donor_nofiltered_1000hvg4.Rdata")
X4 = X ; Y4 = Y %>% as.matrix() ; #Y4 = as.matrix(log1p(Y4)) #; X4 = as(X4, "dgCMatrix") ; Y4 = as(Y4, "dgCMatrix")

load("~/research_multitask/dataset/clr_data_donor_nofiltered_1000hvg5.Rdata")
X5 = X ; Y5 = Y %>% as.matrix(); #Y5 = as.matrix(log1p(Y5)) #; X5 = as(X5, "dgCMatrix") ; Y5 = as(Y5, "dgCMatrix")

load("~/research_multitask/dataset/clr_data_donor_nofiltered_1000hvg6.Rdata")
X6 = X ; Y6 = Y %>% as.matrix() ; #Y6 = as.matrix(log1p(Y6)) #; X6 = as(X6, "dgCMatrix") ; Y6 = as(Y6, "dgCMatrix")

load("~/research_multitask/dataset/clr_data_donor_nofiltered_1000hvg7.Rdata")
X7 = X ; Y7 = Y %>% as.matrix() ; #Y7 = as.matrix(log1p(Y7)) #; X7 = as(X7, "dgCMatrix") ; Y7 = as(Y7, "dgCMatrix")

load("~/research_multitask/dataset/clr_data_donor_nofiltered_1000hvg8.Rdata")
X8 = X ; Y8 = Y %>% as.matrix() ; #Y8 = as.matrix(log1p(Y8)) #; X8 = as(X8, "dgCMatrix") ; Y8 = as(Y8, "dgCMatrix")

rm(X, Y) 

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

p = ncol(X1) ; K = ncol(Y1)

lambda_values = 10^seq(-8,0,length = 10)

train_prop = 0.5 ; val_prop = 0.25 ; test_prop = 0.25

nrep = 20
prediction_mse_vec = rep(NA,nrep)
W_opt_fused_list = list()
tgtidx = 7

#### Model 3) STL B^fused ####

for(j in 1:nrep){
  
  W_opt_fused_list[[j]] = matrix(NA, p, K)
  
  cat('replication = ',j,'\n')
  
  #### this part changes ####
  X0_split_data = split_data(X7, Y7, train_prop, val_prop, test_prop, seed = j-1+nrep*0)
  X1_split_data = split_data(X1, Y1, train_prop, val_prop, test_prop, seed = j-1+nrep*1)
  X2_split_data = split_data(X2, Y2, train_prop, val_prop, test_prop, seed = j-1+nrep*2)
  X3_split_data = split_data(X3, Y3, train_prop, val_prop, test_prop, seed = j-1+nrep*3)
  X4_split_data = split_data(X4, Y4, train_prop, val_prop, test_prop, seed = j-4+nrep*4)
  X5_split_data = split_data(X5, Y5, train_prop, val_prop, test_prop, seed = j-4+nrep*5)
  X6_split_data = split_data(X6, Y6, train_prop, val_prop, test_prop, seed = j-4+nrep*6)
  X7_split_data = split_data(X8, Y8, train_prop, val_prop, test_prop, seed = j-4+nrep*7)
  ####  ####
  
  cv_result = expand.grid(lambda_values)
  cv_result[,2] = NA
  colnames(cv_result) = c('lam1','cvrmse')
  
  X_target = X0_split_data$X_train ; Y_target = X0_split_data$Y_train ; n0 = nrow(X_target)
  X_source1 = X1_split_data$X_train ; Y_source1 = X1_split_data$Y_train ; n1 = nrow(X_source1)
  X_source2 = X2_split_data$X_train ; Y_source2 = X2_split_data$Y_train ; n2 = nrow(X_source2)
  X_source3 = X3_split_data$X_train ; Y_source3 = X3_split_data$Y_train ; n3 = nrow(X_source3)
  X_source4 = X4_split_data$X_train ; Y_source4 = X4_split_data$Y_train ; n4 = nrow(X_source4)
  X_source5 = X5_split_data$X_train ; Y_source5 = X5_split_data$Y_train ; n5 = nrow(X_source5)
  X_source6 = X6_split_data$X_train ; Y_source6 = X6_split_data$Y_train ; n6 = nrow(X_source6)
  X_source7 = X7_split_data$X_train ; Y_source7 = X7_split_data$Y_train ; n7 = nrow(X_source7)
  N = n0+n1+n2+n3+n4+n5+n6+n7
  X0tX0 = eigenMapMatMult(t(X_target), X_target) ; X0tY0 = eigenMapMatMult(t(X_target), Y_target)
  X1tX1 = eigenMapMatMult(t(X_source1), X_source1) ; X1tY1 = eigenMapMatMult(t(X_source1), Y_source1)
  X2tX2 = eigenMapMatMult(t(X_source2), X_source2) ; X2tY2 = eigenMapMatMult(t(X_source2), Y_source2)
  X3tX3 = eigenMapMatMult(t(X_source3), X_source3) ; X3tY3 = eigenMapMatMult(t(X_source3), Y_source3)
  X4tX4 = eigenMapMatMult(t(X_source4), X_source4) ; X4tY4 = eigenMapMatMult(t(X_source4), Y_source4)
  X5tX5 = eigenMapMatMult(t(X_source5), X_source5) ; X5tY5 = eigenMapMatMult(t(X_source5), Y_source5)
  X6tX6 = eigenMapMatMult(t(X_source6), X_source6) ; X6tY6 = eigenMapMatMult(t(X_source6), Y_source6)
  X7tX7 = eigenMapMatMult(t(X_source7), X_source7) ; X7tY7 = eigenMapMatMult(t(X_source7), Y_source7)
  
  
  X_target_ = rbind(X0_split_data$X_train, X0_split_data$X_val) ; Y_target_ = rbind(X0_split_data$Y_train, X0_split_data$Y_val) ; n0_ = nrow(X_target_)
  X_source1_ = rbind(X1_split_data$X_train, X1_split_data$X_val) ; Y_source1_ = rbind(X1_split_data$Y_train, X1_split_data$Y_val) ; n1_ = nrow(X_source1_)
  X_source2_ = rbind(X2_split_data$X_train, X2_split_data$X_val) ; Y_source2_ = rbind(X2_split_data$Y_train, X2_split_data$Y_val) ; n2_ = nrow(X_source2_)
  X_source3_ = rbind(X3_split_data$X_train, X3_split_data$X_val) ; Y_source3_ = rbind(X3_split_data$Y_train, X3_split_data$Y_val) ; n3_ = nrow(X_source3_)
  X_source4_ = rbind(X4_split_data$X_train, X4_split_data$X_val) ; Y_source4_ = rbind(X4_split_data$Y_train, X4_split_data$Y_val) ; n4_ = nrow(X_source4_)
  X_source5_ = rbind(X5_split_data$X_train, X5_split_data$X_val) ; Y_source5_ = rbind(X5_split_data$Y_train, X5_split_data$Y_val) ; n5_ = nrow(X_source5_)
  X_source6_ = rbind(X6_split_data$X_train, X6_split_data$X_val) ; Y_source6_ = rbind(X6_split_data$Y_train, X6_split_data$Y_val) ; n6_ = nrow(X_source6_)
  X_source7_ = rbind(X7_split_data$X_train, X7_split_data$X_val) ; Y_source7_ = rbind(X7_split_data$Y_train, X7_split_data$Y_val) ; n7_ = nrow(X_source7_)
  N_ = n0_+n1_+n2_+n3_+n4_+n5_+n6_+n7_
  X0tX0_ = eigenMapMatMult(t(X_target_), X_target_) ; X0tY0_ = eigenMapMatMult(t(X_target_), Y_target_)
  X1tX1_ = eigenMapMatMult(t(X_source1_), X_source1_) ; X1tY1_ = eigenMapMatMult(t(X_source1_), Y_source1_)
  X2tX2_ = eigenMapMatMult(t(X_source2_), X_source2_) ; X2tY2_ = eigenMapMatMult(t(X_source2_), Y_source2_)
  X3tX3_ = eigenMapMatMult(t(X_source3_), X_source3_) ; X3tY3_ = eigenMapMatMult(t(X_source3_), Y_source3_)
  X4tX4_ = eigenMapMatMult(t(X_source4_), X_source4_) ; X4tY4_ = eigenMapMatMult(t(X_source4_), Y_source4_)
  X5tX5_ = eigenMapMatMult(t(X_source5_), X_source5_) ; X5tY5_ = eigenMapMatMult(t(X_source5_), Y_source5_)
  X6tX6_ = eigenMapMatMult(t(X_source6_), X_source6_) ; X6tY6_ = eigenMapMatMult(t(X_source6_), Y_source6_)
  X7tX7_ = eigenMapMatMult(t(X_source7_), X_source7_) ; X7tY7_ = eigenMapMatMult(t(X_source7_), Y_source7_)
  # X0tX0_ = crossprod(X_target_, X_target_) ; X0tY0_ = crossprod(X_target_,Y_target_)  
  # X1tX1_ = crossprod(X_source1_, X_source1_) ; X1tY1_ = crossprod(X_source1_,Y_source1_)  
  # X2tX2_ = crossprod(X_source2_, X_source2_) ; X2tY2_ = crossprod(X_source2_,Y_source2_)  
  # X3tX3_ = crossprod(X_source3_, X_source3_) ; X3tY3_ = crossprod(X_source3_,Y_source3_)  
  # X4tX4_ = crossprod(X_source4_, X_source4_) ; X4tY4_ = crossprod(X_source4_,Y_source4_)  
  # X5tX5_ = crossprod(X_source5_, X_source5_) ; X5tY5_ = crossprod(X_source5_,Y_source5_)  
  # X6tX6_ = crossprod(X_source6_, X_source6_) ; X6tY6_ = crossprod(X_source6_,Y_source6_)  
  # X7tX7_ = crossprod(X_source7_, X_source7_) ; X7tY7_ = crossprod(X_source7_,Y_source7_) 
  
  
  
  singletask_target_rmse_bytask_transfer = rep(NA, K)
  
  print('Starting task iteration')
  
  for(k in 1:K){ # task loop
    
    if(k%%5==0){cat('k = ',k,'\n')}
    
    X0ty0 = X0tY0[,k] ; X1ty1 = X1tY1[,k] ; X2ty2 = X2tY2[,k] ; X3ty3 = X3tY3[,k] ; X4ty4 = X4tY4[,k] ; X5ty5 = X5tY5[,k] ; X6ty6 = X6tY6[,k] ; X7ty7 = X7tY7[,k]
    # X0ty0 = crossprod(X_source0,Y_source0[,k]) ; X1ty1 = crossprod(X_source1,Y_source1[,k]) ; X2ty2 = crossprod(X_source2,Y_source2[,k]) ; X3ty3 = crossprod(X_source3,Y_source3[,k]) ; X4ty4 = crossprod(X_source4,Y_source4[,k]) ; X5ty5 = crossprod(X_source5,Y_source5[,k]) ; X6ty6 = crossprod(X_source6,Y_source6[,k]) ; X7ty7 = crossprod(X_source7,Y_source7[,k])
    
    cv_result = expand.grid(lambda_values)
    cv_result[,2] = NA
    colnames(cv_result) = c('lam1')
    
    for(i in 1:nrow(cv_result)){ # cv loop
      
      lam1 = cv_result[i,'lam1']
      
      admmm = admm_stl_sevensource(X0 = X_target, y0 = Y_target[,k], X1 = X_source1, y1 = Y_source1[,k], X2 = X_source2, y2 = Y_source2[,k], X3 = X_source3, y3 = Y_source3[,k], 
                                   X4 = X_source4, y4 = Y_source4[,k], X5 = X_source5, y5 = Y_source5[,k], X6 = X_source6, y6 = Y_source6[,k], X7 = X_source7, y7 = Y_source7[,k], 
                                   X0tX0, X0ty0, X1tX1, X1ty1, X2tX2, X2ty2, X3tX3, X3ty3, X4tX4, X4ty4, X5tX5, X5ty5, X6tX6, X6ty6, X7tX7, X7ty7,
                                   lambda0 = lam1, lambda1 = 8*sqrt(n1/N)*lam1, lambda2 = 8*sqrt(n2/N)*lam1, lambda3 = 8*sqrt(n3/N)*lam1, 
                                   lambda4 = 8*sqrt(n4/N)*lam1, lambda5 = 8*sqrt(n5/N)*lam1, lambda6 = 8*sqrt(n6/N)*lam1, lambda7 = 8*sqrt(n7/N)*lam1,
                                   max_iter = 10, tol_prim = 1e-3, tol_dual = 1e-3, verbose = F)
      
      w_opt_fused = n1/N * admmm$Gamma1 + n2/N * admmm$Gamma2 + n3/N * admmm$Gamma3 + n4/N * admmm$Gamma4 + n5/N * admmm$Gamma5 + n6/N * admmm$Gamma6 + n7/N * admmm$Gamma7 + n0/N * admmm$Gamma00 
      
      val_rmse = rmse(X0_split_data$X_val %*% w_opt_fused, X0_split_data$Y_val[,k])
      cv_result[i,'cvrmse'] = val_rmse
      
    } # end cv loop
    
    
    
    
    #### do with best hyperparameter.
    lambda1 = cv_result[which.min(cv_result$cvrmse), 'lam1'] 
    if(k %% 5 == 0){cat('iteration = ',j, 'task = ', k, 'lambda1 = ', lambda1,'\n')}
    
    
    X0ty0_ = X0tY0_[,k] ; X1ty1_ = X1tY1_[,k] ; X2ty2_ = X2tY2_[,k] ; X3ty3_ = X3tY3_[,k] ; X4ty4_ = X4tY4_[,k] ; X5ty5_ = X5tY5_[,k] ; X6ty6_ = X6tY6_[,k] ; X7ty7_ = X7tY7_[,k]
    
    
    admmm = admm_stl_sevensource(X0 = X_target_, y0 = Y_target_[,k], X1 = X_source1_, y1 = Y_source1_[,k], X2 = X_source2_, y2 = Y_source2_[,k], X3 = X_source3_, y3 = Y_source3_[,k], 
                                 X4 = X_source4_, y4 = Y_source4_[,k], X5 = X_source5_, y5 = Y_source5_[,k], X6 = X_source6_, y6 = Y_source6_[,k], X7 = X_source7_, y7 = Y_source7_[,k], 
                                 X0tX0_, X0ty0_, X1tX1_, X1ty1_, X2tX2_, X2ty2_, X3tX3_, X3ty3_, X4tX4_, X4ty4_, X5tX5_, X5ty5_, X6tX6_, X6ty6_, X7tX7_, X7ty7_,
                                 lambda0 = lambda1, lambda1 = 8*sqrt(n1_/N_)*lambda1, lambda2 = 8*sqrt(n2_/N_)*lambda1, lambda3 = 8*sqrt(n3_/N_)*lambda1, 
                                 lambda4 = 8*sqrt(n4_/N_)*lambda1, lambda5 = 8*sqrt(n5_/N_)*lambda1, lambda6 = 8*sqrt(n6_/N_)*lambda1, lambda7 = 8*sqrt(n7_/N_)*lambda1,
                                 max_iter = 50, tol_prim = 1e-3, tol_dual = 1e-3, verbose = F)
    
    w_opt_fused = n1_/N_ * admmm$Gamma1 + n2_/N_ * admmm$Gamma2 + n3_/N_ * admmm$Gamma3 + n4_/N_ * admmm$Gamma4 + n5_/N_ * admmm$Gamma5 + n6_/N_ * admmm$Gamma6 + n7_/N_ * admmm$Gamma7 + n0_/N_ * admmm$Gamma00 
    test_rmse_k = rmse(X0_split_data$X_test %*% w_opt_fused, X0_split_data$Y_test[,k])
    singletask_target_rmse_bytask_transfer[k] = test_rmse_k
    
    W_opt_fused_list[[j]][,k] = w_opt_fused
    
  } # end task loop
  
  prediction_mse_vec[j] = sqrt(mean(singletask_target_rmse_bytask_transfer^2))
  
  save(prediction_mse_vec, W_opt_fused_list, file = paste("~/research_multitask/results/summary_logtransform_target",tgtidx,".Rdata",sep=""))
}# end iteration loop

save(prediction_mse_vec, W_opt_fused_list, file = paste("~/research_multitask/results/summary_logtransform_target",tgtidx,".Rdata",sep=""))

