
library(glmnet)
library(Matrix)
library(dplyr)
source('functions_joint.R')
source('simulationsetting.R')

tgtidx = 1

load(paste('../clr_data_donor_nofiltered_500hvg1.Rdata',sep=""))
X1 = X ; Y1 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg2.Rdata',sep=""))
X2 = X ; Y2 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg3.Rdata',sep=""))
X3 = X ; Y3 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg4.Rdata',sep=""))
X4 = X ; Y4 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg5.Rdata',sep=""))
X5 = X ; Y5 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg6.Rdata',sep=""))
X6 = X ; Y6 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg7.Rdata',sep=""))
X7 = X ; Y7 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg8.Rdata',sep=""))
X8 = X ; Y8 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg9.Rdata',sep=""))
X9 = X ; Y9 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg10.Rdata',sep=""))
X10 = X ; Y10 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg11.Rdata',sep=""))
X11 = X ; Y11 = as.matrix(Y)
load(paste('../clr_data_donor_nofiltered_500hvg12.Rdata',sep=""))
X12 = X ; Y12 = as.matrix(Y)
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
lambda0_vec = rep(NA,nrep)
W_opt_fused_list = list()

for(j in 1:nrep){
  
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
  
  X_target = X0_split_data$X_train ; Y_target = X0_split_data$Y_train
  X_source1 = X1_split_data$X_train ; Y_source1 = X1_split_data$Y_train
  X_source2 = X2_split_data$X_train ; Y_source2 = X2_split_data$Y_train
  X_source3 = X3_split_data$X_train ; Y_source3 = X3_split_data$Y_train
  X_source4 = X4_split_data$X_train ; Y_source4 = X4_split_data$Y_train
  X_source5 = X5_split_data$X_train ; Y_source5 = X5_split_data$Y_train
  X_source6 = X6_split_data$X_train ; Y_source6 = X6_split_data$Y_train
  X_source7 = X7_split_data$X_train ; Y_source7 = X7_split_data$Y_train
  X_source8 = X8_split_data$X_train ; Y_source8 = X8_split_data$Y_train
  X_source9 = X9_split_data$X_train ; Y_source9 = X9_split_data$Y_train
  X_source10 = X10_split_data$X_train ; Y_source10 = X10_split_data$Y_train
  X_source11 = X11_split_data$X_train ; Y_source11 = X11_split_data$Y_train
  
  X0tX0 = crossprod(X_target, X_target) ; X0tY0 = crossprod(X_target,Y_target) ; n0 = nrow(X_target)
  X1tX1 = crossprod(X_source1, X_source1) ; X1tY1 = crossprod(X_source1,Y_source1) ; n1 = nrow(X_source1)
  X2tX2 = crossprod(X_source2, X_source2) ; X2tY2 = crossprod(X_source2,Y_source2) ; n2 = nrow(X_source2)
  X3tX3 = crossprod(X_source3, X_source3) ; X3tY3 = crossprod(X_source3,Y_source3) ; n3 = nrow(X_source3)
  X4tX4 = crossprod(X_source4, X_source4) ; X4tY4 = crossprod(X_source4,Y_source4) ; n4 = nrow(X_source4)
  X5tX5 = crossprod(X_source5, X_source5) ; X5tY5 = crossprod(X_source5,Y_source5) ; n5 = nrow(X_source5)
  X6tX6 = crossprod(X_source6, X_source6) ; X6tY6 = crossprod(X_source6,Y_source6) ; n6 = nrow(X_source6)
  X7tX7 = crossprod(X_source7, X_source7) ; X7tY7 = crossprod(X_source7,Y_source7) ; n7 = nrow(X_source7)
  X8tX8 = crossprod(X_source8, X_source8) ; X8tY8 = crossprod(X_source8,Y_source8) ; n8 = nrow(X_source8)
  X9tX9 = crossprod(X_source9, X_source9) ; X9tY9 = crossprod(X_source9,Y_source9) ; n9 = nrow(X_source9)
  X10tX10 = crossprod(X_source10, X_source10) ; X10tY10 = crossprod(X_source10,Y_source10) ; n10 = nrow(X_source10)
  X11tX11 = crossprod(X_source11, X_source11) ; X11tY11 = crossprod(X_source11,Y_source11) ; n11 = nrow(X_source11)
  N = n0+n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11
  
  for(i in 1:nrow(cv_result)){ # start cv
    
    cat('lambda = ', lambda_values[i], '\n')
    
    lam1 = cv_result[i,'lam1']
    
    admmm = admm_eleven_source(X0 = X_target, Y0 = Y_target, X1 = X_source1, Y1 = Y_source1, X2 = X_source2, Y2 = Y_source2, X3 = X_source3, Y3 = Y_source3, X4 = X_source4, Y4 = Y_source4, X5 = X_source5, Y5 = Y_source5, X6 = X_source6, Y6 = Y_source6, X7 = X_source7, Y7 = Y_source7, X8 = X_source8, Y8 = Y_source8, X9 = X_source9, Y9 = Y_source9, X10 = X_source10, Y10 = Y_source10, X11 = X_source11, Y11 = Y_source11,
                               X0tX0, X0tY0, X1tX1, X1tY1, X2tX2, X2tY2, X3tX3, X3tY3, X4tX4, X4tY4, X5tX5, X5tY5, X6tX6, X6tY6, X7tX7, X7tY7, X8tX8, X8tY8, X9tX9, X9tY9, X10tX10, X10tY10, X11tX11, X11tY11,
                               lambda0 = lam1, lambda1 = 8*sqrt(n1/N)*lam1, lambda2 = 8*sqrt(n2/N)*lam1, lambda3 = 8*sqrt(n3/N)*lam1, lambda4 = 8*sqrt(n4/N)*lam1, lambda5 = 8*sqrt(n5/N)*lam1, lambda6 = 8*sqrt(n6/N)*lam1, lambda7 = 8*sqrt(n7/N)*lam1, lambda8 = 8*sqrt(n8/N)*lam1, lambda9 = 8*sqrt(n9/N)*lam1, lambda10 = 8*sqrt(n10/N)*lam1, lambda11 = 8*sqrt(n11/N)*lam1,
                               max_iter = max_iter, tol_prim = tol_prim, tol_dual = tol_dual, verbose = F)
    
    W_opt_fused = n1/N * admmm$Gamma1 + n2/N * admmm$Gamma2 + n3/N * admmm$Gamma3 + n4/N * admmm$Gamma4 + n5/N * admmm$Gamma5 + n6/N * admmm$Gamma6 + n7/N * admmm$Gamma7 + n8/N * admmm$Gamma8 + n9/N * admmm$Gamma9 + n10/N * admmm$Gamma10 + n11/N * admmm$Gamma11 + n0/N * admmm$Gamma00
    
    val_rmse = rmse(X0_split_data$X_val %*% W_opt_fused, X0_split_data$Y_val)
    cv_result[i,'cvrmse'] = val_rmse
    
  }
  
  print(cv_result)
  
  lambda1 = cv_result[which.min(cv_result$cvrmse), 'lam1']
  cat('iteration = ',j, 'lambda1 = ', lambda1,'\n')
  
  X_target = rbind(X0_split_data$X_train, X0_split_data$X_val) ; Y_target = rbind(X0_split_data$Y_train, X0_split_data$Y_val) ; n0 = nrow(X_target)
  X_source1 = rbind(X1_split_data$X_train, X1_split_data$X_val) ; Y_source1 = rbind(X1_split_data$Y_train, X1_split_data$Y_val) ; n1 = nrow(X_source1)
  X_source2 = rbind(X2_split_data$X_train, X2_split_data$X_val) ; Y_source2 = rbind(X2_split_data$Y_train, X2_split_data$Y_val) ; n2 = nrow(X_source2)
  X_source3 = rbind(X3_split_data$X_train, X3_split_data$X_val) ; Y_source3 = rbind(X3_split_data$Y_train, X3_split_data$Y_val) ; n3 = nrow(X_source3)
  X_source4 = rbind(X4_split_data$X_train, X4_split_data$X_val) ; Y_source4 = rbind(X4_split_data$Y_train, X4_split_data$Y_val) ; n4 = nrow(X_source4)
  X_source5 = rbind(X5_split_data$X_train, X5_split_data$X_val) ; Y_source5 = rbind(X5_split_data$Y_train, X5_split_data$Y_val) ; n5 = nrow(X_source5)
  X_source6 = rbind(X6_split_data$X_train, X6_split_data$X_val) ; Y_source6 = rbind(X6_split_data$Y_train, X6_split_data$Y_val) ; n6 = nrow(X_source6)
  X_source7 = rbind(X7_split_data$X_train, X7_split_data$X_val) ; Y_source7 = rbind(X7_split_data$Y_train, X7_split_data$Y_val) ; n7 = nrow(X_source7)
  X_source8 = rbind(X8_split_data$X_train, X8_split_data$X_val) ; Y_source8 = rbind(X8_split_data$Y_train, X8_split_data$Y_val) ; n8 = nrow(X_source8)
  X_source9 = rbind(X9_split_data$X_train, X9_split_data$X_val) ; Y_source9 = rbind(X9_split_data$Y_train, X9_split_data$Y_val) ; n9 = nrow(X_source9)
  X_source10 = rbind(X10_split_data$X_train, X10_split_data$X_val) ; Y_source10 = rbind(X10_split_data$Y_train, X10_split_data$Y_val) ; n10 = nrow(X_source10)
  X_source11 = rbind(X11_split_data$X_train, X11_split_data$X_val) ; Y_source11 = rbind(X11_split_data$Y_train, X11_split_data$Y_val) ; n11 = nrow(X_source11)
  
  X0tX0 = crossprod(X_target, X_target) ; X0tY0 = crossprod(X_target,Y_target) ; n0 = nrow(X_target)
  X1tX1 = crossprod(X_source1, X_source1) ; X1tY1 = crossprod(X_source1,Y_source1) ; n1 = nrow(X_source1)
  X2tX2 = crossprod(X_source2, X_source2) ; X2tY2 = crossprod(X_source2,Y_source2) ; n2 = nrow(X_source2)
  X3tX3 = crossprod(X_source3, X_source3) ; X3tY3 = crossprod(X_source3,Y_source3) ; n3 = nrow(X_source3)
  X4tX4 = crossprod(X_source4, X_source4) ; X4tY4 = crossprod(X_source4,Y_source4) ; n4 = nrow(X_source4)
  X5tX5 = crossprod(X_source5, X_source5) ; X5tY5 = crossprod(X_source5,Y_source5) ; n5 = nrow(X_source5)
  X6tX6 = crossprod(X_source6, X_source6) ; X6tY6 = crossprod(X_source6,Y_source6) ; n6 = nrow(X_source6)
  X7tX7 = crossprod(X_source7, X_source7) ; X7tY7 = crossprod(X_source7,Y_source7) ; n7 = nrow(X_source7)
  X8tX8 = crossprod(X_source8, X_source8) ; X8tY8 = crossprod(X_source8,Y_source8) ; n8 = nrow(X_source8)
  X9tX9 = crossprod(X_source9, X_source9) ; X9tY9 = crossprod(X_source9,Y_source9) ; n9 = nrow(X_source9)
  X10tX10 = crossprod(X_source10, X_source10) ; X10tY10 = crossprod(X_source10,Y_source10) ; n10 = nrow(X_source10)
  X11tX11 = crossprod(X_source11, X_source11) ; X11tY11 = crossprod(X_source11,Y_source11) ; n11 = nrow(X_source11)
  N = n0+n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11
  
  admmm = admm_eleven_source(X0 = X_target, Y0 = Y_target, X1 = X_source1, Y1 = Y_source1, X2 = X_source2, Y2 = Y_source2, X3 = X_source3, Y3 = Y_source3, X4 = X_source4, Y4 = Y_source4, X5 = X_source5, Y5 = Y_source5, X6 = X_source6, Y6 = Y_source6, X7 = X_source7, Y7 = Y_source7, X8 = X_source8, Y8 = Y_source8, X9 = X_source9, Y9 = Y_source9, X10 = X_source10, Y10 = Y_source10, X11 = X_source11, Y11 = Y_source11,
                             X0tX0, X0tY0, X1tX1, X1tY1, X2tX2, X2tY2, X3tX3, X3tY3, X4tX4, X4tY4, X5tX5, X5tY5, X6tX6, X6tY6, X7tX7, X7tY7, X8tX8, X8tY8, X9tX9, X9tY9, X10tX10, X10tY10, X11tX11, X11tY11,
                             lambda0 = lambda1, lambda1 = 8*sqrt(n1/N)*lambda1, lambda2 = 8*sqrt(n2/N)*lambda1, lambda3 = 8*sqrt(n3/N)*lambda1, lambda4 = 8*sqrt(n4/N)*lambda1, lambda5 = 8*sqrt(n5/N)*lambda1, lambda6 = 8*sqrt(n6/N)*lambda1, lambda7 = 8*sqrt(n7/N)*lambda1, lambda8 = 8*sqrt(n8/N)*lambda1, lambda9 = 8*sqrt(n9/N)*lambda1, lambda10 = 8*sqrt(n10/N)*lambda1, lambda11 = 8*sqrt(n11/N)*lambda1,
                             max_iter = max_iter, tol_prim = tol_prim, tol_dual = tol_dual, verbose = T)
  
  cat('iter =', admmm$iter, '/', max_iter, ' rprim =', tail(na.omit(admmm$r_primal00),1), ' sdual =', tail(na.omit(admmm$s_dual00),1), '\n')
  
  W_opt_fused = n1/N * admmm$Gamma1 + n2/N * admmm$Gamma2 + n3/N * admmm$Gamma3 + n4/N * admmm$Gamma4 + n5/N * admmm$Gamma5 + n6/N * admmm$Gamma6 + n7/N * admmm$Gamma7 + n8/N * admmm$Gamma8 + n9/N * admmm$Gamma9 + n10/N * admmm$Gamma10 + n11/N * admmm$Gamma11 + n0/N * admmm$Gamma00
  
  parameter_list = W_opt_fused
  
  test_predictions = sweep(as.matrix(X0_split_data$X_test %*% W_opt_fused), 2, X0_split_data$muY, "+")
  Y_test_true      = sweep(as.matrix(X0_split_data$Y_test), 2, X0_split_data$muY, "+")
  
  lambda0_vec[j] = lambda1
  prediction_mse_vec[j] = sse(test_predictions, Y_test_true)
  W_opt_fused_list[[j]] = W_opt_fused
  
  print(lambda0_vec)
  print(prediction_mse_vec)
  cat(100*mean(rowMeans(W_opt_fused) == 0),'% of the features are zero','\n')
  
  save(prediction_mse_vec, lambda0_vec, W_opt_fused_list, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
}
save(prediction_mse_vec, lambda0_vec, W_opt_fused_list, file = paste("summary_logtransform_target",tgtidx,".Rdata",sep=""))
