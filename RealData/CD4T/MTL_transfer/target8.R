library(glmnet)
library(Matrix)
library(dplyr)

set.seed(1)

source('~/Library/CloudStorage/GoogleDrive-96limtotoro@gmail.com/My Drive/research_multitask/functions.R')

load(paste('../clr_data_donor_nofiltered_1000hvg1.Rdata',sep=""))
X1 = X ; Y1 = as.matrix(Y) ; #Y1 = as.matrix(log1p(Y1)) #; X1 = as(X1, "dgCMatrix") ; Y1 = as(Y1, "dgCMatrix")

load(paste('../clr_data_donor_nofiltered_1000hvg2.Rdata',sep=""))
X2 = X ; Y2 = as.matrix(Y)  ; #Y2 = as.matrix(log1p(Y2)) #; X2 = as(X2, "dgCMatrix") ; Y2 = as(Y2, "dgCMatrix")

load(paste('../clr_data_donor_nofiltered_1000hvg3.Rdata',sep=""))
X3 = X ; Y3 = as.matrix(Y)  ; #Y3 = as.matrix(log1p(Y3)) #; X3 = as(X3, "dgCMatrix") ; Y3 = as(Y3, "dgCMatrix")

load(paste('../clr_data_donor_nofiltered_1000hvg4.Rdata',sep=""))
X4 = X ; Y4 = as.matrix(Y)  ; #Y4 = as.matrix(log1p(Y4)) #; X4 = as(X4, "dgCMatrix") ; Y4 = as(Y4, "dgCMatrix")

load(paste('../clr_data_donor_nofiltered_1000hvg5.Rdata',sep=""))
X5 = X ; Y5 = as.matrix(Y)  ; #Y5 = as.matrix(log1p(Y5)) #; X5 = as(X5, "dgCMatrix") ; Y5 = as(Y5, "dgCMatrix")

load(paste('../clr_data_donor_nofiltered_1000hvg6.Rdata',sep=""))
X6 = X ; Y6 = as.matrix(Y)  ; #Y6 = as.matrix(log1p(Y6)) #; X6 = as(X6, "dgCMatrix") ; Y6 = as(Y6, "dgCMatrix")

load(paste('../clr_data_donor_nofiltered_1000hvg7.Rdata',sep=""))
X7 = X ; Y7 = as.matrix(Y)  ; #Y7 = as.matrix(log1p(Y7)) #; X7 = as(X7, "dgCMatrix") ; Y7 = as(Y7, "dgCMatrix")

load(paste('../clr_data_donor_nofiltered_1000hvg8.Rdata',sep=""))
X8 = X ; Y8 = as.matrix(Y)  ; #Y8 = as.matrix(log1p(Y8)) #; X8 = as(X8, "dgCMatrix") ; Y8 = as(Y8, "dgCMatrix")

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


prediction_mse_vec = rep(NA,nrep)
lambda0_vec = rep(NA,nrep)
W_opt_fused_list = list()

for(j in 1:nrep){
  
  cat('replication = ',j,'\n')
  
  #### this part changes ####
  X0_split_data = split_data(X8, Y8, train_prop, val_prop, test_prop, seed = j-1+nrep*0)
  X1_split_data = split_data(X1, Y1, train_prop, val_prop, test_prop, seed = j-1+nrep*1)
  X2_split_data = split_data(X2, Y2, train_prop, val_prop, test_prop, seed = j-1+nrep*2)
  X3_split_data = split_data(X3, Y3, train_prop, val_prop, test_prop, seed = j-1+nrep*3)
  X4_split_data = split_data(X4, Y4, train_prop, val_prop, test_prop, seed = j-4+nrep*4)
  X5_split_data = split_data(X5, Y5, train_prop, val_prop, test_prop, seed = j-4+nrep*5)
  X6_split_data = split_data(X6, Y6, train_prop, val_prop, test_prop, seed = j-4+nrep*6)
  X7_split_data = split_data(X7, Y7, train_prop, val_prop, test_prop, seed = j-4+nrep*7)
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
  
  X0tX0 = crossprod(X_target, X_target) ; X0tY0 = crossprod(X_target,Y_target) ; n0 = nrow(X_target)
  X1tX1 = crossprod(X_source1, X_source1) ; X1tY1 = crossprod(X_source1,Y_source1) ; n1 = nrow(X_source1)
  X2tX2 = crossprod(X_source2, X_source2) ; X2tY2 = crossprod(X_source2,Y_source2) ; n2 = nrow(X_source2)
  X3tX3 = crossprod(X_source3, X_source3) ; X3tY3 = crossprod(X_source3,Y_source3) ; n3 = nrow(X_source3)
  X4tX4 = crossprod(X_source4, X_source4) ; X4tY4 = crossprod(X_source4,Y_source4) ; n4 = nrow(X_source4)
  X5tX5 = crossprod(X_source5, X_source5) ; X5tY5 = crossprod(X_source5,Y_source5) ; n5 = nrow(X_source5)
  X6tX6 = crossprod(X_source6, X_source6) ; X6tY6 = crossprod(X_source6,Y_source6) ; n6 = nrow(X_source6)
  X7tX7 = crossprod(X_source7, X_source7) ; X7tY7 = crossprod(X_source7,Y_source7) ; n7 = nrow(X_source7)
  N = n0+n1+n2+n3+n4+n5+n6+n7
  
  for(i in 1:nrow(cv_result)){ # start cv
    
    cat('lambda = ', lambda_values[i], '\n')
    
    lam1 = cv_result[i,'lam1'] #; lam2 = cv_result[i,'lam2'] ; lam3 = cv_result[i,'lam3']
    
    admmm = admm_seven_source(X0 = X_target, Y0 = Y_target, X1 = X_source1, Y1 = Y_source1, X2 = X_source2, Y2 = Y_source2, X3 = X_source3, Y3 = Y_source3, X4 = X_source4, Y4 = Y_source4, X5 = X_source5, Y5 = Y_source5, X6 = X_source6, Y6 = Y_source6, X7 = X_source7, Y7 = Y_source7,
                              X0tX0, X0tY0, X1tX1, X1tY1, X2tX2, X2tY2, X3tX3, X3tY3, X4tX4, X4tY4, X5tX5, X5tY5, X6tX6, X6tY6, X7tX7, X7tY7,
                              lambda0 = lam1, lambda1 = 8*sqrt(n1/N)*lam1, lambda2 = 8*sqrt(n2/N)*lam1, lambda3 = 8*sqrt(n3/N)*lam1, lambda4 = 8*sqrt(n4/N)*lam1, lambda5 = 8*sqrt(n5/N)*lam1, lambda6 = 8*sqrt(n6/N)*lam1, lambda7 = 8*sqrt(n7/N)*lam1,
                              max_iter = 400, tol_prim = sqrt(K)*5e-5, tol_dual = sqrt(K)*5e-5, verbose = F)
    
    W_opt_fused = n1/N * admmm$Gamma1 + n2/N * admmm$Gamma2 + n3/N * admmm$Gamma3 + n4/N * admmm$Gamma4 + n5/N * admmm$Gamma5 + n6/N * admmm$Gamma6 + n7/N * admmm$Gamma7 + n0/N * admmm$Gamma00
    
    val_rmse = rmse(X0_split_data$X_val %*% W_opt_fused, X0_split_data$Y_val)
    cv_result[i,'cvrmse'] = val_rmse
    
  }
  
  lambda1 = cv_result[which.min(cv_result$cvrmse), 'lam1'] ; # lambda2 = cv_result[which.min(cv_result$cvrmse), 'lam2'] ; lambda3 = cv_result[which.min(cv_result$cvrmse), 'lam3']
  cat('iteration = ',j, 'lambda1 = ', lambda1,'\n')#, ', lambda2 = ', lambda2, 'lambda3 = ', lambda3, '\n')
  
  X_target = rbind(X0_split_data$X_train, X0_split_data$X_val) ; Y_target = rbind(X0_split_data$Y_train, X0_split_data$Y_val) ; n0 = nrow(X_target)
  X_source1 = rbind(X1_split_data$X_train, X1_split_data$X_val) ; Y_source1 = rbind(X1_split_data$Y_train, X1_split_data$Y_val) ; n1 = nrow(X_source1)
  X_source2 = rbind(X2_split_data$X_train, X2_split_data$X_val) ; Y_source2 = rbind(X2_split_data$Y_train, X2_split_data$Y_val) ; n2 = nrow(X_source2)
  X_source3 = rbind(X3_split_data$X_train, X3_split_data$X_val) ; Y_source3 = rbind(X3_split_data$Y_train, X3_split_data$Y_val) ; n3 = nrow(X_source3)
  X_source4 = rbind(X4_split_data$X_train, X4_split_data$X_val) ; Y_source4 = rbind(X4_split_data$Y_train, X4_split_data$Y_val) ; n4 = nrow(X_source4)
  X_source5 = rbind(X5_split_data$X_train, X5_split_data$X_val) ; Y_source5 = rbind(X5_split_data$Y_train, X5_split_data$Y_val) ; n5 = nrow(X_source5)
  X_source6 = rbind(X6_split_data$X_train, X6_split_data$X_val) ; Y_source6 = rbind(X6_split_data$Y_train, X6_split_data$Y_val) ; n6 = nrow(X_source6)
  X_source7 = rbind(X7_split_data$X_train, X7_split_data$X_val) ; Y_source7 = rbind(X7_split_data$Y_train, X7_split_data$Y_val) ; n7 = nrow(X_source7)
  
  X0tX0 = crossprod(X_target, X_target) ; X0tY0 = crossprod(X_target,Y_target) ; n0 = nrow(X_target)
  X1tX1 = crossprod(X_source1, X_source1) ; X1tY1 = crossprod(X_source1,Y_source1) ; n1 = nrow(X_source1)
  X2tX2 = crossprod(X_source2, X_source2) ; X2tY2 = crossprod(X_source2,Y_source2) ; n2 = nrow(X_source2)
  X3tX3 = crossprod(X_source3, X_source3) ; X3tY3 = crossprod(X_source3,Y_source3) ; n3 = nrow(X_source3)
  X4tX4 = crossprod(X_source4, X_source4) ; X4tY4 = crossprod(X_source4,Y_source4) ; n4 = nrow(X_source4)
  X5tX5 = crossprod(X_source5, X_source5) ; X5tY5 = crossprod(X_source5,Y_source5) ; n5 = nrow(X_source5)
  X6tX6 = crossprod(X_source6, X_source6) ; X6tY6 = crossprod(X_source6,Y_source6) ; n6 = nrow(X_source6)
  X7tX7 = crossprod(X_source7, X_source7) ; X7tY7 = crossprod(X_source7,Y_source7) ; n7 = nrow(X_source7)
  N = n0+n1+n2+n3+n4+n5+n6+n7
  
  admmm = admm_seven_source(X0 = X_target, Y0 = Y_target, X1 = X_source1, Y1 = Y_source1, X2 = X_source2, Y2 = Y_source2, X3 = X_source3, Y3 = Y_source3, X4 = X_source4, Y4 = Y_source4, X5 = X_source5, Y5 = Y_source5, X6 = X_source6, Y6 = Y_source6, X7 = X_source7, Y7 = Y_source7,
                            X0tX0, X0tY0, X1tX1, X1tY1, X2tX2, X2tY2, X3tX3, X3tY3, X4tX4, X4tY4, X5tX5, X5tY5, X6tX6, X6tY6, X7tX7, X7tY7,
                            lambda0 = lambda1, lambda1 = 8*sqrt(n1/N)*lambda1, lambda2 = 8*sqrt(n2/N)*lambda1, lambda3 = 8*sqrt(n3/N)*lambda1, lambda4 = 8*sqrt(n4/N)*lambda1, lambda5 = 8*sqrt(n5/N)*lambda1, lambda6 = 8*sqrt(n6/N)*lambda1, lambda7 = 8*sqrt(n7/N)*lambda1,
                            max_iter = 400, tol_prim = sqrt(K)*5e-5, tol_dual = sqrt(K)*5e-5, verbose = T)
  
  W_opt_fused = n1/N * admmm$Gamma1 + n2/N * admmm$Gamma2 + n3/N * admmm$Gamma3 + n4/N * admmm$Gamma4 + n5/N * admmm$Gamma5 + n6/N * admmm$Gamma6 + n7/N * admmm$Gamma7 + n0/N * admmm$Gamma00
  
  parameter_list = W_opt_fused
  
  test_predictions = X0_split_data$X_test %*% W_opt_fused ; Y_test_true = X0_split_data$Y_test
  
  
  lambda0_vec[j] = lambda1
  prediction_mse_vec[j] = sse(test_predictions, X0_split_data$Y_test)
  W_opt_fused_list[[j]] = W_opt_fused
  
  print(lambda0_vec)
  print(prediction_mse_vec)
  cat(100*mean(rowMeans(W_opt_fused) == 0),'% of the features are zero','\n')
  save(prediction_mse_vec, lambda0_vec, W_opt_fused_list, file = "summary_logtransform_target8.Rdata")
}

save(prediction_mse_vec, lambda0_vec, W_opt_fused_list, file = "summary_logtransform_target8.Rdata")