library(glmnet)
library(Matrix)
source('simulationsetting.R')
library(dplyr)
tgtidx = 1
source('functions_stl_cluster.R')
source('simulationsetting.R')

load('../clr_data_donor_nofiltered_500hvg1.Rdata')
X1 = X ; Y1 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg2.Rdata')
X2 = X ; Y2 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg3.Rdata')
X3 = X ; Y3 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg4.Rdata')
X4 = X ; Y4 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg5.Rdata')
X5 = X ; Y5 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg6.Rdata')
X6 = X ; Y6 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg7.Rdata')
X7 = X ; Y7 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg8.Rdata')
X8 = X ; Y8 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg9.Rdata')
X9 = X ; Y9 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg10.Rdata')
X10 = X ; Y10 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg11.Rdata')
X11 = X ; Y11 = Y %>% as.matrix()
load('../clr_data_donor_nofiltered_500hvg12.Rdata')
X12 = X ; Y12 = Y %>% as.matrix()
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

p = ncol(X1) ; K = ncol(Y1)
prediction_mse_vec = rep(NA,nrep)
W_opt_fused_list = list()

#### Model 3) STL B^fused ####
for(j in 1:nrep){
  
  W_opt_fused_list[[j]] = matrix(NA, p, K)
  
  set.seed(j)
  
  cat('replication = ',j,'\n')
  
  #### this part changes ####
  X0_split_data = center_XY(split_data(X1, Y1, train_prop, val_prop, test_prop))
  X1_split_data = center_XY(split_data(X2, Y2, train_prop, val_prop, test_prop))
  X2_split_data = center_XY(split_data(X3, Y3, train_prop, val_prop, test_prop))
  X3_split_data = center_XY(split_data(X4, Y4, train_prop, val_prop, test_prop))
  X4_split_data = center_XY(split_data(X5, Y5, train_prop, val_prop, test_prop))
  X5_split_data = center_XY(split_data(X6, Y6, train_prop, val_prop, test_prop))
  X6_split_data = center_XY(split_data(X7, Y7, train_prop, val_prop, test_prop))
  X7_split_data = center_XY(split_data(X8, Y8, train_prop, val_prop, test_prop))
  X8_split_data = center_XY(split_data(X9, Y9, train_prop, val_prop, test_prop))
  X9_split_data = center_XY(split_data(X10, Y10, train_prop, val_prop, test_prop))
  X10_split_data = center_XY(split_data(X11, Y11, train_prop, val_prop, test_prop))
  X11_split_data = center_XY(split_data(X12, Y12, train_prop, val_prop, test_prop))
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
  X_source8 = X8_split_data$X_train ; Y_source8 = X8_split_data$Y_train ; n8 = nrow(X_source8)
  X_source9 = X9_split_data$X_train ; Y_source9 = X9_split_data$Y_train ; n9 = nrow(X_source9)
  X_source10 = X10_split_data$X_train ; Y_source10 = X10_split_data$Y_train ; n10 = nrow(X_source10)
  X_source11 = X11_split_data$X_train ; Y_source11 = X11_split_data$Y_train ; n11 = nrow(X_source11)
  N = n0+n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11
  X0tX0 = eigenMapMatMult(t(X_target), X_target) ; X0tY0 = eigenMapMatMult(t(X_target), Y_target)
  X1tX1 = eigenMapMatMult(t(X_source1), X_source1) ; X1tY1 = eigenMapMatMult(t(X_source1), Y_source1)
  X2tX2 = eigenMapMatMult(t(X_source2), X_source2) ; X2tY2 = eigenMapMatMult(t(X_source2), Y_source2)
  X3tX3 = eigenMapMatMult(t(X_source3), X_source3) ; X3tY3 = eigenMapMatMult(t(X_source3), Y_source3)
  X4tX4 = eigenMapMatMult(t(X_source4), X_source4) ; X4tY4 = eigenMapMatMult(t(X_source4), Y_source4)
  X5tX5 = eigenMapMatMult(t(X_source5), X_source5) ; X5tY5 = eigenMapMatMult(t(X_source5), Y_source5)
  X6tX6 = eigenMapMatMult(t(X_source6), X_source6) ; X6tY6 = eigenMapMatMult(t(X_source6), Y_source6)
  X7tX7 = eigenMapMatMult(t(X_source7), X_source7) ; X7tY7 = eigenMapMatMult(t(X_source7), Y_source7)
  X8tX8 = eigenMapMatMult(t(X_source8), X_source8) ; X8tY8 = eigenMapMatMult(t(X_source8), Y_source8)
  X9tX9 = eigenMapMatMult(t(X_source9), X_source9) ; X9tY9 = eigenMapMatMult(t(X_source9), Y_source9)
  X10tX10 = eigenMapMatMult(t(X_source10), X_source10) ; X10tY10 = eigenMapMatMult(t(X_source10), Y_source10)
  X11tX11 = eigenMapMatMult(t(X_source11), X_source11) ; X11tY11 = eigenMapMatMult(t(X_source11), Y_source11)
  
  
  X_target_ = rbind(X0_split_data$X_train, X0_split_data$X_val) ; Y_target_ = rbind(X0_split_data$Y_train, X0_split_data$Y_val) ; n0_ = nrow(X_target_)
  X_source1_ = rbind(X1_split_data$X_train, X1_split_data$X_val) ; Y_source1_ = rbind(X1_split_data$Y_train, X1_split_data$Y_val) ; n1_ = nrow(X_source1_)
  X_source2_ = rbind(X2_split_data$X_train, X2_split_data$X_val) ; Y_source2_ = rbind(X2_split_data$Y_train, X2_split_data$Y_val) ; n2_ = nrow(X_source2_)
  X_source3_ = rbind(X3_split_data$X_train, X3_split_data$X_val) ; Y_source3_ = rbind(X3_split_data$Y_train, X3_split_data$Y_val) ; n3_ = nrow(X_source3_)
  X_source4_ = rbind(X4_split_data$X_train, X4_split_data$X_val) ; Y_source4_ = rbind(X4_split_data$Y_train, X4_split_data$Y_val) ; n4_ = nrow(X_source4_)
  X_source5_ = rbind(X5_split_data$X_train, X5_split_data$X_val) ; Y_source5_ = rbind(X5_split_data$Y_train, X5_split_data$Y_val) ; n5_ = nrow(X_source5_)
  X_source6_ = rbind(X6_split_data$X_train, X6_split_data$X_val) ; Y_source6_ = rbind(X6_split_data$Y_train, X6_split_data$Y_val) ; n6_ = nrow(X_source6_)
  X_source7_ = rbind(X7_split_data$X_train, X7_split_data$X_val) ; Y_source7_ = rbind(X7_split_data$Y_train, X7_split_data$Y_val) ; n7_ = nrow(X_source7_)
  X_source8_ = rbind(X8_split_data$X_train, X8_split_data$X_val) ; Y_source8_ = rbind(X8_split_data$Y_train, X8_split_data$Y_val) ; n8_ = nrow(X_source8_)
  X_source9_ = rbind(X9_split_data$X_train, X9_split_data$X_val) ; Y_source9_ = rbind(X9_split_data$Y_train, X9_split_data$Y_val) ; n9_ = nrow(X_source9_)
  X_source10_ = rbind(X10_split_data$X_train, X10_split_data$X_val) ; Y_source10_ = rbind(X10_split_data$Y_train, X10_split_data$Y_val) ; n10_ = nrow(X_source10_)
  X_source11_ = rbind(X11_split_data$X_train, X11_split_data$X_val) ; Y_source11_ = rbind(X11_split_data$Y_train, X11_split_data$Y_val) ; n11_ = nrow(X_source11_)
  N_ = n0_+n1_+n2_+n3_+n4_+n5_+n6_+n7_+n8_+n9_+n10_+n11_
  X0tX0_ = eigenMapMatMult(t(X_target_), X_target_) ; X0tY0_ = eigenMapMatMult(t(X_target_), Y_target_)
  X1tX1_ = eigenMapMatMult(t(X_source1_), X_source1_) ; X1tY1_ = eigenMapMatMult(t(X_source1_), Y_source1_)
  X2tX2_ = eigenMapMatMult(t(X_source2_), X_source2_) ; X2tY2_ = eigenMapMatMult(t(X_source2_), Y_source2_)
  X3tX3_ = eigenMapMatMult(t(X_source3_), X_source3_) ; X3tY3_ = eigenMapMatMult(t(X_source3_), Y_source3_)
  X4tX4_ = eigenMapMatMult(t(X_source4_), X_source4_) ; X4tY4_ = eigenMapMatMult(t(X_source4_), Y_source4_)
  X5tX5_ = eigenMapMatMult(t(X_source5_), X_source5_) ; X5tY5_ = eigenMapMatMult(t(X_source5_), Y_source5_)
  X6tX6_ = eigenMapMatMult(t(X_source6_), X_source6_) ; X6tY6_ = eigenMapMatMult(t(X_source6_), Y_source6_)
  X7tX7_ = eigenMapMatMult(t(X_source7_), X_source7_) ; X7tY7_ = eigenMapMatMult(t(X_source7_), Y_source7_)
  X8tX8_ = eigenMapMatMult(t(X_source8_), X_source8_) ; X8tY8_ = eigenMapMatMult(t(X_source8_), Y_source8_)
  X9tX9_ = eigenMapMatMult(t(X_source9_), X_source9_) ; X9tY9_ = eigenMapMatMult(t(X_source9_), Y_source9_)
  X10tX10_ = eigenMapMatMult(t(X_source10_), X_source10_) ; X10tY10_ = eigenMapMatMult(t(X_source10_), Y_source10_)
  X11tX11_ = eigenMapMatMult(t(X_source11_), X_source11_) ; X11tY11_ = eigenMapMatMult(t(X_source11_), Y_source11_)
  
  
  singletask_target_rmse_bytask_transfer = rep(NA, K)
  
  print('Starting task iteration')
  
  for(k in 1:K){ # task loop
    
    if(k%%5==0){cat('k = ',k,'\n')}
    
    X0ty0 = X0tY0[,k] ; X1ty1 = X1tY1[,k] ; X2ty2 = X2tY2[,k] ; X3ty3 = X3tY3[,k] ; X4ty4 = X4tY4[,k] ; X5ty5 = X5tY5[,k] ; X6ty6 = X6tY6[,k] ; X7ty7 = X7tY7[,k] ; X8ty8 = X8tY8[,k] ; X9ty9 = X9tY9[,k] ; X10ty10 = X10tY10[,k] ; X11ty11 = X11tY11[,k]
    cv_result = expand.grid(lambda_values)
    cv_result[,2] = NA
    colnames(cv_result) = c('lam1','cvrmse')
    
    for(i in 1:nrow(cv_result)){ # cv loop
      
      lam1 = cv_result[i,'lam1']
      
      admmm = admm_stl_elevensource(X0 = X_target, y0 = Y_target[,k], X1 = X_source1, y1 = Y_source1[,k], X2 = X_source2, y2 = Y_source2[,k], X3 = X_source3, y3 = Y_source3[,k], X4 = X_source4, y4 = Y_source4[,k], X5 = X_source5, y5 = Y_source5[,k], X6 = X_source6, y6 = Y_source6[,k], X7 = X_source7, y7 = Y_source7[,k], X8 = X_source8, y8 = Y_source8[,k], X9 = X_source9, y9 = Y_source9[,k], X10 = X_source10, y10 = Y_source10[,k], X11 = X_source11, y11 = Y_source11[,k],
                                    X0tX0, X0ty0, X1tX1, X1ty1, X2tX2, X2ty2, X3tX3, X3ty3, X4tX4, X4ty4, X5tX5, X5ty5, X6tX6, X6ty6, X7tX7, X7ty7, X8tX8, X8ty8, X9tX9, X9ty9, X10tX10, X10ty10, X11tX11, X11ty11,
                                    lambda0 = lam1, lambda1 = 8*sqrt(n1/N)*lam1, lambda2 = 8*sqrt(n2/N)*lam1, lambda3 = 8*sqrt(n3/N)*lam1, lambda4 = 8*sqrt(n4/N)*lam1, lambda5 = 8*sqrt(n5/N)*lam1, lambda6 = 8*sqrt(n6/N)*lam1, lambda7 = 8*sqrt(n7/N)*lam1, lambda8 = 8*sqrt(n8/N)*lam1, lambda9 = 8*sqrt(n9/N)*lam1, lambda10 = 8*sqrt(n10/N)*lam1, lambda11 = 8*sqrt(n11/N)*lam1,
                                    max_iter = max_iter_stl, tol_prim = tol_prim_stl, tol_dual = tol_dual_stl, verbose = F)
      
      w_opt_fused = n1/N * admmm$Gamma1 + n2/N * admmm$Gamma2 + n3/N * admmm$Gamma3 + n4/N * admmm$Gamma4 + n5/N * admmm$Gamma5 + n6/N * admmm$Gamma6 + n7/N * admmm$Gamma7 + n8/N * admmm$Gamma8 + n9/N * admmm$Gamma9 + n10/N * admmm$Gamma10 + n11/N * admmm$Gamma11 + n0/N * admmm$Gamma00
      
      val_rmse = rmse(X0_split_data$X_val %*% w_opt_fused, X0_split_data$Y_val[,k])
      cv_result[i,'cvrmse'] = val_rmse
      
    } # end cv loop
    
    #### do with best hyperparameter.
    lambda1 = cv_result[which.min(cv_result$cvrmse), 'lam1']
    if(k %% 5 == 0){cat('iteration = ',j, 'task = ', k, 'lambda1 = ', lambda1,'\n')}
    
    X0ty0_ = X0tY0_[,k] ; X1ty1_ = X1tY1_[,k] ; X2ty2_ = X2tY2_[,k] ; X3ty3_ = X3tY3_[,k] ; X4ty4_ = X4tY4_[,k] ; X5ty5_ = X5tY5_[,k] ; X6ty6_ = X6tY6_[,k] ; X7ty7_ = X7tY7_[,k] ; X8ty8_ = X8tY8_[,k] ; X9ty9_ = X9tY9_[,k] ; X10ty10_ = X10tY10_[,k] ; X11ty11_ = X11tY11_[,k]
    
    admmm = admm_stl_elevensource(X0 = X_target_, y0 = Y_target_[,k], X1 = X_source1_, y1 = Y_source1_[,k], X2 = X_source2_, y2 = Y_source2_[,k], X3 = X_source3_, y3 = Y_source3_[,k], X4 = X_source4_, y4 = Y_source4_[,k], X5 = X_source5_, y5 = Y_source5_[,k], X6 = X_source6_, y6 = Y_source6_[,k], X7 = X_source7_, y7 = Y_source7_[,k], X8 = X_source8_, y8 = Y_source8_[,k], X9 = X_source9_, y9 = Y_source9_[,k], X10 = X_source10_, y10 = Y_source10_[,k], X11 = X_source11_, y11 = Y_source11_[,k],
                                  X0tX0_, X0ty0_, X1tX1_, X1ty1_, X2tX2_, X2ty2_, X3tX3_, X3ty3_, X4tX4_, X4ty4_, X5tX5_, X5ty5_, X6tX6_, X6ty6_, X7tX7_, X7ty7_, X8tX8_, X8ty8_, X9tX9_, X9ty9_, X10tX10_, X10ty10_, X11tX11_, X11ty11_,
                                  lambda0 = lambda1, lambda1 = 8*sqrt(n1_/N_)*lambda1, lambda2 = 8*sqrt(n2_/N_)*lambda1, lambda3 = 8*sqrt(n3_/N_)*lambda1, lambda4 = 8*sqrt(n4_/N_)*lambda1, lambda5 = 8*sqrt(n5_/N_)*lambda1, lambda6 = 8*sqrt(n6_/N_)*lambda1, lambda7 = 8*sqrt(n7_/N_)*lambda1, lambda8 = 8*sqrt(n8_/N_)*lambda1, lambda9 = 8*sqrt(n9_/N_)*lambda1, lambda10 = 8*sqrt(n10_/N_)*lambda1, lambda11 = 8*sqrt(n11_/N_)*lambda1,
                                  max_iter = max_iter_stl, tol_prim = tol_prim_stl, tol_dual = tol_dual_stl, verbose = F)
    
    if(k %% 10 == 0){cat('  k=',k,' iter=',admmm$iter,'/',max_iter_stl,'\n')}
    
    w_opt_fused = n1_/N_ * admmm$Gamma1 + n2_/N_ * admmm$Gamma2 + n3_/N_ * admmm$Gamma3 + n4_/N_ * admmm$Gamma4 + n5_/N_ * admmm$Gamma5 + n6_/N_ * admmm$Gamma6 + n7_/N_ * admmm$Gamma7 + n8_/N_ * admmm$Gamma8 + n9_/N_ * admmm$Gamma9 + n10_/N_ * admmm$Gamma10 + n11_/N_ * admmm$Gamma11 + n0_/N_ * admmm$Gamma00
    test_rmse_k = rmse(as.numeric(X0_split_data$X_test %*% w_opt_fused) + X0_split_data$muY[k],
                       X0_split_data$Y_test[,k] + X0_split_data$muY[k])
    singletask_target_rmse_bytask_transfer[k] = test_rmse_k
    
    W_opt_fused_list[[j]][,k] = w_opt_fused
    
  } # end task loop
  
  prediction_mse_vec[j] = sqrt(mean(singletask_target_rmse_bytask_transfer^2))
  
  save(prediction_mse_vec, W_opt_fused_list, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
}# end iteration loop
save(prediction_mse_vec, W_opt_fused_list, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
