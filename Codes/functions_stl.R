library(CVXR)
library(glmnet)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('~/Library/CloudStorage/GoogleDrive-96limtotoro@gmail.com/My Drive/research_multitask/MatMult_mac.cpp')

split_data <- function(X, Y, train_prop = 0.6, val_prop = 0.2, test_prop = 0.2) {
  #set.seed(seed)
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

rmse = function(x,y){
  sqrt(mean((x-y)^2))
}
sse = function(x,y){
  sum((x-y)^2)
}
# oracle_translasso <- function(X0, y0, X1, y1, 
#                               lambda0, lambda1){
#   
#   n_target = nrow(X0)
#   n_source = nrow(X1)
#   
#   model = glmnet(X1, y1, alpha = 1, 
#                  lambda = lambda1, family = "gaussian") 
#   w1 = as.numeric(model$beta)
#   
#   model = glmnet(X0, y0 - X0 %*% w1, alpha = 1,
#                  lambda = lambda0, family = "gaussian")
#   delta = as.numeric(model$beta)
#   
#   beta = w1 + delta
#   
#   return(beta)
# }
# 
# oracle_translasso_twosource <- function(X0, y0, X1, y1, X2, y2,
#                                         lambda0, lambda1){
#                                           
#   n_target = nrow(X0) 
#   n_source = nrow(X1) + nrow(X2)
#   
#   X_source = rbind(X1,X2)
#   y_source = c(y1,y2)
#   
#   model = glmnet(X_source, y_source, alpha = 1, 
#                  lambda = lambda1, family = "gaussian") 
#   w1 = as.numeric(model$beta)
#   
#   model = glmnet(X0, y0 - X0 %*% w1, alpha = 1,
#                  lambda = lambda0, family = "gaussian")
#   delta = as.numeric(model$beta)
#   
#   beta = w1 + delta
#   
#   return(beta)
# }
# 
# oracle_translasso_threesource <- function(X0, y0, X1, y1, X2, y2, X3, y3, 
#                                           lambda0, lambda1){
#                               
#   
#   n_target = nrow(X0) 
#   n_source = nrow(X1) + nrow(X2) + nrow(X3)
#   #p <- ncol(X0)
#   
#   X_source = rbind(X1,X2,X3)
#   y_source = c(y1,y2,y3)
#   
#   model = glmnet(X_source, y_source, alpha = 1, 
#                  lambda = lambda1, family = "gaussian") 
#   w1 = as.numeric(model$beta)
#   
#   model = glmnet(X0, y0 - X0 %*% w1, alpha = 1,
#                  lambda = lambda0, family = "gaussian")
#   delta = as.numeric(model$beta)
#   
#   beta = w1 + delta
#   
#   return(beta)
# }


soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

# cvx_stl = function(X0,y0,X1,y1,
#                    lambda0, lambda1){ # CVXR implementation of He et al. 2024 method with 1/(2N) weight for all dataset
#   
#   N = nrow(X0) + nrow(X1)
#   p = ncol(X0)
#   
#   w_source = Variable(p,1) ; w_target = Variable(p,1)
#   loss = 1/(2*N)*(sum_squares(X0 %*% w_target - y0) + sum_squares(X1 %*% w_source - y1))
#   objective = Minimize(loss + lambda0 * p_norm(w_target, p = 1) + lambda1 * p_norm(w_target - w_source, p = 1))
#   problem = Problem(objective)
#   result = solve(problem)
#   w_opt_target = as.matrix(result$getValue(w_target))
#   w_opt_source = as.matrix(result$getValue(w_source))
#   
#   return(list(B0 = w_opt_target, B1 = w_opt_source, result = result))
# }


admm_stl <- function(X0, y0, X1, y1, 
                     X0tX0 = NULL, X0ty0 = NULL, X1tX1 = NULL, X1ty1 = NULL,
                     lambda0, lambda1, 
                     rho00 = 1, rho01 = 1, rho1 = 1, 
                     mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                     max_iter = 2000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){ # ADMM implementation of He et al. 2024 method with 1/(2N) weight for all dataset
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1)
  
  objvalue = rep(NA, max_iter)
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho1value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma1 <- matrix(0, p, 1)
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  V1 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0) ; X1tX1 = crossprod(X1,X1)
    X0ty0 = crossprod(X0,y0) ; X1ty1 = crossprod(X1,y1)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma1_old = Gamma1
    
    #### Step 1: B0 update ####
    lhs = 1/N*X0tX0 + (rho00 + rho01)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1)
    B0 <- solve(lhs, rhs)
    
    #### Step 2: B1 update ####
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### Step 3: Gamma00 update ####
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Step 4: Gamma01 update ####
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    
    #### Step 5: Gamma1 update ####
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)

    #### Step 6–8: Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    V1 <- V1 + B1 - Gamma1
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)

    if (verbose && iter %% 500 == 0) {
      cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f \n",
                  iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal1,base=10), log(s_dual1,base=10) ))
    }
    
    objvalue[iter] = 1/(2*N)*(l2norm_sq(y0 - X0 %*% Gamma00) + l2norm_sq(y1 - X1 %*% Gamma1)) + lambda0 * sum(abs(Gamma00)) + lambda1 * sum(abs(Gamma1-Gamma00)) # objective value
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho1value[iter] = rho1
    
    if (max(r_primal00, r_primal01, r_primal1) < tol_prim && max(s_dual00, s_dual01, s_dual1) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, 
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma1 = Gamma1,
              U0 = U0, U1 = U1, V1 = V1, 
              objvalue = objvalue, 
              rho00 = rho00value, rho01 = rho01value, rho1 = rho1value,
              r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal1 = r_primal1value, s_dual1 = s_dual1value))
}

admm_stl_twosource <- function(X0, y0, X1, y1, X2, y2,
                               X0tX0 = NULL, X0ty0 = NULL, X1tX1 = NULL, X1ty1 = NULL, X2tX2 = NULL, X2ty2 = NULL,
                               lambda0, lambda1, lambda2,
                               rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1,
                               mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                               max_iter = 2000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){ # ADMM implementation of He et al. 2024 method with 1/(2N) weight for all dataset
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1) + nrow(X2)
  
  objvalue = rep(NA, max_iter)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho02value = rep(NA, max_iter)
  
  rho1value = rep(NA, max_iter)
  rho2value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  B2 <- matrix(0, p, 1)
  
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma02 <- matrix(0, p, 1)
  
  Gamma1 <- matrix(0, p, 1)
  Gamma2 <- matrix(0, p, 1)
  
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  U2 <- matrix(0, p, 1)
  
  V1 <- matrix(0, p, 1)
  V2 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0) ; X1tX1 = crossprod(X1,X1) ; X2tX2 = crossprod(X2,X2) 
    X0ty0 = crossprod(X0,y0) ; X1ty1 = crossprod(X1,y1) ; X2ty2 = crossprod(X2,y2) 
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    
    #### Step 1: B0 update ####
    lhs = 1/N*X0tX0 + (rho00 + rho01 + rho02)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2)
    B0 <- solve(lhs, rhs)
    
    #### Step 2: B1 update ####
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### Step 3: B2 update ####
    lhs = 1/N*X2tX2 + rho2*Ip
    rhs = 1/N*X2ty2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### Step 4: Gamma00 update ####
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Step 5: Gamma01 update ####
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    
    #### Step 6: Gamma02 update ####
    Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
    
    #### Step 7: Gamma1 update ####
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
    
    #### Step 8: Gamma2 update ####
    Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
    
    #### Step 9–13: Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
    
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
    r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
    
    if (verbose && iter %% 500 == 0) {
      cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid02: %.2f | Dual Resid02: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f | Primal Resid2: %.2f | Dual Resid2: %.2f \n",
                  iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10)))
    }
    objvalue[iter] = 1/(2*N)*(l2norm_sq(y0 - X0 %*% Gamma00) + l2norm_sq(y1 - X1 %*% Gamma1) + l2norm_sq(y2 - X2 %*% Gamma2)) + lambda0 * sum(abs(Gamma00)) + lambda1 * sum(abs(Gamma1-Gamma00)) + lambda2 * sum(abs(Gamma2-Gamma00)) # objective value
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
    
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal02 > mu * s_dual02){
      rho02 = min(rho02 * tau_incr, rho_max)
      U2 = U2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual02 > mu * r_primal02){
      rho02 = rho02 / tau_decr
      U2 = U2 * tau_decr # change of variable required for the scaled variable
    }
    
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal2 > mu * s_dual2){
      rho2 = min(rho2 * tau_incr, rho_max)
      V2 = V2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual2 > mu * r_primal2){
      rho2 = rho2 / tau_decr
      V2 = V2 * tau_decr # change of variable required for the scaled variable
    }
    
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho02value[iter] = rho02
    
    rho1value[iter] = rho1
    rho2value[iter] = rho2
    
    if (max(r_primal00, r_primal01, r_primal02, r_primal1, r_primal2) < tol_prim && max(s_dual00, s_dual01, s_dual02, s_dual1, s_dual2) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, B2 = B2,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02,
              Gamma1 = Gamma1, Gamma2 = Gamma2,
              U0 = U0, U1 = U1, U2 = U2,
              V1 = V1, V2 = V2,
              objvalue = objvalue, 
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value,
              rho1 = rho1value, rho2 = rho2value,
              
              r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value))
}

admm_stl_threesource <- function(X0, y0, X1, y1, X2, y2, X3, y3,
                                 X0tX0 = NULL, X0ty0 = NULL, X1tX1 = NULL, X1ty1 = NULL, X2tX2 = NULL, X2ty2 = NULL, X3tX3 = NULL, X3ty3 = NULL,
                                 lambda0, lambda1, lambda2, lambda3,
                                 rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho1 = 1, rho2 = 1, rho3 = 1,
                                 mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                 max_iter = 5000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){ # ADMM implementation of He et al. 2024 method with 1/(2N) weight for all dataset
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3)
  
  objvalue = rep(NA, max_iter)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho02value = rep(NA, max_iter)
  rho03value = rep(NA, max_iter)
  
  rho1value = rep(NA, max_iter)
  rho2value = rep(NA, max_iter)
  rho3value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
  
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  B2 <- matrix(0, p, 1)
  B3 <- matrix(0, p, 1)
  
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma02 <- matrix(0, p, 1)
  Gamma03 <- matrix(0, p, 1)
  
  Gamma1 <- matrix(0, p, 1)
  Gamma2 <- matrix(0, p, 1)
  Gamma3 <- matrix(0, p, 1)
  
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  U2 <- matrix(0, p, 1)
  U3 <- matrix(0, p, 1)
  
  V1 <- matrix(0, p, 1)
  V2 <- matrix(0, p, 1)
  V3 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0) ; X1tX1 = crossprod(X1,X1) ; X2tX2 = crossprod(X2,X2) ; X3tX3 = crossprod(X3,X3)
    X0ty0 = crossprod(X0,y0) ; X1ty1 = crossprod(X1,y1) ; X2ty2 = crossprod(X2,y2) ; X3ty3 = crossprod(X3,y3)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    # if(iter %% 100 == 0){
    #   cat('iter = ', iter, '\n')
    # }
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma03_old = Gamma03
    
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    Gamma3_old = Gamma3
    
    #### Step 1: B0 update ####
    #cat('step 1 \n')
    lhs = 1/N*X0tX0 + (rho00 + rho01 + rho02 + rho03)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3)
    B0 <- solve(lhs, rhs)
    #print(dim(B0))
    
    #### Step 2: B1 update ####
    #cat('step 2 \n')
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    #print(dim(B1))
    
    #### Step 3: B2 update ####
    #cat('step 3 \n')
    lhs = 1/N*X2tX2 + rho2*Ip
    rhs = 1/N*X2ty2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    #print(dim(B2))
    
    #### Step 4: B3 update ####
    #cat('step 4 \n')
    lhs = 1/N*X3tX3 + rho3*Ip
    rhs = 1/N*X3ty3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    #print(dim(B3))
    
    #### Step 5: Gamma00 update ####
    #cat('step 5 \n')
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Step 6: Gamma01 update ####
    #cat('step 6 \n')
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    
    #### Step 7: Gamma02 update ####
    #cat('step 7 \n')
    Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
    
    #### Step 8: Gamma03 update ####
    #cat('step 8 \n')
    Gamma03 = Gamma3 + soft_threshold(B0+U3-Gamma3, lambda3/rho03)
    
    #### Step 9: Gamma1 update ####
    #cat('step 9 \n')
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
    
    #### Step 10: Gamma2 update ####
    #cat('step 10 \n')
    Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
    
    #### Step 11: Gamma3 update ####
    #cat('step 11 \n')
    Gamma3 = Gamma03 + soft_threshold(B3+V3-Gamma03, lambda3/rho3)
    
    #### Step 12–18: Dual updates ####
    #cat('step 12\n')
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    U3 <- U3 + B0 - Gamma03
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    V3 <- V3 + B3 - Gamma3
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
    r_primal03 = l2norm(B0 - Gamma03) ; s_dual03 = rho03 * l2norm(Gamma03 - Gamma03_old)
    
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
    r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
    r_primal3 = l2norm(B3 - Gamma3) ; s_dual3 = rho3 * l2norm(Gamma3 - Gamma3_old)
    
    if (verbose && iter %% 1000 == 0) {
      cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid02: %.2f | Dual Resid02: %.2f | Primal Resid03: %.2f |Dual Resid03: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f | Primal Resid2: %.2f | Dual Resid2: %.2f | Primal Resid3: %.2f | Dual Resid3: %.2f \n",
                  iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal03,base=10), log(s_dual03,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10), log(r_primal3,base=10), log(s_dual3,base=10)))
    }
    
    objvalue[iter] = 1/(2*N)*(l2norm_sq(y0 - X0 %*% B0) + l2norm_sq(y1 - X1 %*% B1) + l2norm_sq(y2 - X2 %*% B2) + l2norm_sq(y3 - X3 %*% B3)) + lambda0 * sum(abs(Gamma00)) + lambda1 * sum(abs(Gamma1-Gamma00)) + lambda2 * sum(abs(Gamma2-Gamma00)) + lambda3 * sum(abs(Gamma3-Gamma00)) # objective value
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
    r_primal03value[iter] = r_primal03 ; s_dual03value[iter] = s_dual03
    
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
    r_primal3value[iter] = r_primal3 ; s_dual3value[iter] = s_dual3
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal02 > mu * s_dual02){
      rho02 = min(rho02 * tau_incr, rho_max)
      U2 = U2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual02 > mu * r_primal02){
      rho02 = rho02 / tau_decr
      U2 = U2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal03 > mu * s_dual03){
      rho03 = min(rho03 * tau_incr, rho_max)
      U3 = U3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual03 > mu * r_primal03){
      rho03 = rho03 / tau_decr
      U3 = U3 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal2 > mu * s_dual2){
      rho2 = min(rho2 * tau_incr, rho_max)
      V2 = V2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual2 > mu * r_primal2){
      rho2 = rho2 / tau_decr
      V2 = V2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal3 > mu * s_dual3){
      rho3 = min(rho3 * tau_incr, rho_max)
      V3 = V3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual3 > mu * r_primal3){
      rho3 = rho3 / tau_decr
      V3 = V3 * tau_decr # change of variable required for the scaled variable
    }
    
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho02value[iter] = rho02
    rho03value[iter] = rho03
    
    rho1value[iter] = rho1
    rho2value[iter] = rho2
    rho3value[iter] = rho3
    
    if (max(r_primal00, r_primal01, r_primal02, r_primal03, r_primal1, r_primal2, r_primal3) < tol_prim && max(s_dual00, s_dual01, s_dual02, s_dual03, s_dual1, s_dual2, s_dual3) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, 
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, 
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, 
              V1 = V1, V2 = V2, V3 = V3, 
              objvalue = objvalue, 
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, 
              rho1 = rho1value, rho2 = rho2value, rho3 = rho3value,
              
              r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value))
}

admm_stl_fivesource <- function(X0, y0, X1, y1, X2, y2, X3, y3, X4, y4, X5, y5,
                                X0tX0 = NULL, X0ty0 = NULL, X1tX1 = NULL, X1ty1 = NULL, X2tX2 = NULL, X2ty2 = NULL, X3tX3 = NULL, X3ty3 = NULL, X4tX4 = NULL, X4ty4 = NULL, X5tX5 = NULL, X5ty5 = NULL,
                                lambda0, lambda1, lambda2, lambda3,lambda4,lambda5,
                                rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1,
                                mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                max_iter = 5000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){ # ADMM implementation of He et al. 2024 method with 1/(2N) weight for all dataset
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3) + nrow(X4) + nrow(X5)
  
  objvalue = rep(NA, max_iter)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho02value = rep(NA, max_iter)
  rho03value = rep(NA, max_iter)
  rho04value = rep(NA, max_iter)
  rho05value = rep(NA, max_iter)
  
  rho1value = rep(NA, max_iter)
  rho2value = rep(NA, max_iter)
  rho3value = rep(NA, max_iter)
  rho4value = rep(NA, max_iter)
  rho5value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
  r_primal04value = rep(NA, max_iter) ; s_dual04value = rep(NA, max_iter)
  r_primal05value = rep(NA, max_iter) ; s_dual05value = rep(NA, max_iter)
  
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
  r_primal4value = rep(NA, max_iter) ; s_dual4value = rep(NA, max_iter)
  r_primal5value = rep(NA, max_iter) ; s_dual5value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  B2 <- matrix(0, p, 1)
  B3 <- matrix(0, p, 1)
  B4 <- matrix(0, p, 1)
  B5 <- matrix(0, p, 1)
  
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma02 <- matrix(0, p, 1)
  Gamma03 <- matrix(0, p, 1)
  Gamma04 <- matrix(0, p, 1)
  Gamma05 <- matrix(0, p, 1)
  
  Gamma1 <- matrix(0, p, 1)
  Gamma2 <- matrix(0, p, 1)
  Gamma3 <- matrix(0, p, 1)
  Gamma4 <- matrix(0, p, 1)
  Gamma5 <- matrix(0, p, 1)
  
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  U2 <- matrix(0, p, 1)
  U3 <- matrix(0, p, 1)
  U4 <- matrix(0, p, 1)
  U5 <- matrix(0, p, 1)
  
  V1 <- matrix(0, p, 1)
  V2 <- matrix(0, p, 1)
  V3 <- matrix(0, p, 1)
  V4 <- matrix(0, p, 1)
  V5 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0) ; X1tX1 = crossprod(X1,X1) ; X2tX2 = crossprod(X2,X2) ; X3tX3 = crossprod(X3,X3) ; X4tX4 = crossprod(X4,X4) ; X5tX5 = crossprod(X5,X5)
    X0ty0 = crossprod(X0,y0) ; X1ty1 = crossprod(X1,y1) ; X2ty2 = crossprod(X2,y2) ; X3ty3 = crossprod(X3,y3) ; X4ty4 = crossprod(X4,y4) ; X5ty5 = crossprod(X5,y5)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    # if(iter %% 100 == 0){
    #   cat('iter = ', iter, '\n')
    # }
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma03_old = Gamma03
    Gamma04_old = Gamma04
    Gamma05_old = Gamma05
    
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    Gamma3_old = Gamma3
    Gamma4_old = Gamma4
    Gamma5_old = Gamma5
    
    #### Step 1: B0 update ####
    lhs = 1/N*X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5)
    B0 <- solve(lhs, rhs)
    
    #### Step 2: B1 update ####
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### Step 3: B2 update ####
    lhs = 1/N*X2tX2 + rho2*Ip
    rhs = 1/N*X2ty2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### Step 4: B3 update ####
    lhs = 1/N*X3tX3 + rho3*Ip
    rhs = 1/N*X3ty3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### Step 5: B4 update ####
    lhs = 1/N*X4tX4 + rho4*Ip
    rhs = 1/N*X4ty4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### Step 6: B5 update ####
    lhs = 1/N*X5tX5 + rho5*Ip
    rhs = 1/N*X5ty5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    
    #### Step 7: Gamma00 update ####
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Step 8: Gamma01 update ####
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    
    #### Step 9: Gamma02 update ####
    Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
    
    #### Step 10: Gamma03 update ####
    Gamma03 = Gamma3 + soft_threshold(B0+U3-Gamma3, lambda3/rho03)
    
    #### Step 11: Gamma04 update ####
    Gamma04 = Gamma4 + soft_threshold(B0+U4-Gamma4, lambda4/rho04)
    
    #### Step 12: Gamma05 update ####
    Gamma05 = Gamma5 + soft_threshold(B0+U5-Gamma5, lambda5/rho05)
    
    
    
    #### Step 13: Gamma1 update ####
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
    
    #### Step 14: Gamma2 update ####
    Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
    
    #### Step 15: Gamma3 update ####
    Gamma3 = Gamma03 + soft_threshold(B3+V3-Gamma03, lambda3/rho3)
    
    #### Step 16: Gamma4 update ####
    Gamma4 = Gamma04 + soft_threshold(B4+V4-Gamma04, lambda4/rho4)
    
    #### Step 17: Gamma5 update ####
    Gamma5 = Gamma05 + soft_threshold(B5+V5-Gamma05, lambda5/rho5)
    
    
    #### Step 18–28: Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    U3 <- U3 + B0 - Gamma03
    U4 <- U4 + B0 - Gamma04
    U5 <- U5 + B0 - Gamma05
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    V3 <- V3 + B3 - Gamma3
    V4 <- V4 + B4 - Gamma4
    V5 <- V5 + B5 - Gamma5
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
    r_primal03 = l2norm(B0 - Gamma03) ; s_dual03 = rho03 * l2norm(Gamma03 - Gamma03_old)
    r_primal04 = l2norm(B0 - Gamma04) ; s_dual04 = rho04 * l2norm(Gamma04 - Gamma04_old)
    r_primal05 = l2norm(B0 - Gamma05) ; s_dual05 = rho05 * l2norm(Gamma05 - Gamma05_old)
    
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
    r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
    r_primal3 = l2norm(B3 - Gamma3) ; s_dual3 = rho3 * l2norm(Gamma3 - Gamma3_old)
    r_primal4 = l2norm(B4 - Gamma4) ; s_dual4 = rho4 * l2norm(Gamma4 - Gamma4_old)
    r_primal5 = l2norm(B5 - Gamma5) ; s_dual5 = rho5 * l2norm(Gamma5 - Gamma5_old)
    
    if (verbose && iter %% 1000 == 0) {
      cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid02: %.2f | Dual Resid02: %.2f | Primal Resid03: %.2f |Dual Resid03: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f | Primal Resid2: %.2f | Dual Resid2: %.2f | Primal Resid3: %.2f | Dual Resid3: %.2f \n",
                  iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal03,base=10), log(s_dual03,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10), log(r_primal3,base=10), log(s_dual3,base=10)))
    }
    
    objvalue[iter] = 1/(2*N)*(l2norm_sq(y0 - X0 %*% B0) + l2norm_sq(y1 - X1 %*% B1) + l2norm_sq(y2 - X2 %*% B2) + l2norm_sq(y3 - X3 %*% B3) + l2norm_sq(y4 - X4 %*% B4) + l2norm_sq(y5 - X5 %*% B5)) +
      lambda0 * sum(abs(Gamma00)) + lambda1 * sum(abs(Gamma1-Gamma00)) + lambda2 * sum(abs(Gamma2-Gamma00)) + lambda3 * sum(abs(Gamma3-Gamma00)) + lambda4 * sum(abs(Gamma4-Gamma00)) + lambda5 * sum(abs(Gamma5-Gamma00))  
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
    r_primal03value[iter] = r_primal03 ; s_dual03value[iter] = s_dual03
    r_primal04value[iter] = r_primal04 ; s_dual04value[iter] = s_dual04
    r_primal05value[iter] = r_primal05 ; s_dual05value[iter] = s_dual05
    
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
    r_primal3value[iter] = r_primal3 ; s_dual3value[iter] = s_dual3
    r_primal4value[iter] = r_primal4 ; s_dual4value[iter] = s_dual4
    r_primal5value[iter] = r_primal5 ; s_dual5value[iter] = s_dual5
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal02 > mu * s_dual02){
      rho02 = min(rho02 * tau_incr, rho_max)
      U2 = U2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual02 > mu * r_primal02){
      rho02 = rho02 / tau_decr
      U2 = U2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal03 > mu * s_dual03){
      rho03 = min(rho03 * tau_incr, rho_max)
      U3 = U3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual03 > mu * r_primal03){
      rho03 = rho03 / tau_decr
      U3 = U3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal04 > mu * s_dual04){
      rho04 = min(rho04 * tau_incr, rho_max)
      U4 = U4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual04 > mu * r_primal04){
      rho04 = rho04 / tau_decr
      U4 = U4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal05 > mu * s_dual05){
      rho05 = min(rho05 * tau_incr, rho_max)
      U5 = U5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual05 > mu * r_primal05){
      rho05 = rho05 / tau_decr
      U5 = U5 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal2 > mu * s_dual2){
      rho2 = min(rho2 * tau_incr, rho_max)
      V2 = V2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual2 > mu * r_primal2){
      rho2 = rho2 / tau_decr
      V2 = V2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal3 > mu * s_dual3){
      rho3 = min(rho3 * tau_incr, rho_max)
      V3 = V3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual3 > mu * r_primal3){
      rho3 = rho3 / tau_decr
      V3 = V3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal4 > mu * s_dual4){
      rho4 = min(rho4 * tau_incr, rho_max)
      V4 = V4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual4 > mu * r_primal4){
      rho4 = rho4 / tau_decr
      V4 = V4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal5 > mu * s_dual5){
      rho5 = min(rho5 * tau_incr, rho_max)
      V5 = V5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual5 > mu * r_primal5){
      rho5 = rho5 / tau_decr
      V5 = V5 * tau_decr # change of variable required for the scaled variable
    }
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho02value[iter] = rho02
    rho03value[iter] = rho03
    rho04value[iter] = rho04
    rho05value[iter] = rho05
    
    rho1value[iter] = rho1
    rho2value[iter] = rho2
    rho3value[iter] = rho3
    rho4value[iter] = rho4
    rho5value[iter] = rho5
    
    if (max(r_primal00, r_primal01, r_primal02, r_primal03, r_primal04, r_primal05, r_primal1, r_primal2, r_primal3, r_primal4, r_primal5) < tol_prim && max(s_dual00, s_dual01, s_dual02, s_dual03, s_dual04, s_dual05, s_dual1, s_dual2, s_dual3, s_dual4, s_dual5) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5,
              objvalue = objvalue, 
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, rho04 = rho04value, rho05 = rho05value,
              rho1 = rho1value, rho2 = rho2value, rho3 = rho3value, rho4 = rho4value, rho5 = rho5value,
              
              r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              r_primal04 = r_primal04value, s_dual04 = s_dual04value,
              r_primal05 = r_primal05value, s_dual05 = s_dual05value,
              
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value,
              r_primal4 = r_primal4value, s_dual4 = s_dual4value,
              r_primal5 = r_primal5value, s_dual5 = s_dual5value))
}

admm_stl_sevensource <- function(X0, y0, X1, y1, X2, y2, X3, y3, X4, y4, X5, y5, X6, y6, X7, y7,
                                 X0tX0 = NULL, X0ty0 = NULL, 
                                 X1tX1 = NULL, X1ty1 = NULL, X2tX2 = NULL, X2ty2 = NULL, X3tX3 = NULL, X3ty3 = NULL, X4tX4 = NULL, X4ty4 = NULL, 
                                 X5tX5 = NULL, X5ty5 = NULL, X6tX6 = NULL, X6ty6 = NULL, X7tX7 = NULL, X7ty7 = NULL,
                                 lambda0, lambda1, lambda2, lambda3,lambda4,lambda5,lambda6,lambda7,
                                 rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1, rho07 = 1,
                                 rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1, rho7 = 1,
                                 mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                 max_iter = 5000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3) + nrow(X4) + nrow(X5) + nrow(X6) + nrow(X7)
  
  #objvalue = rep(NA, max_iter)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho02value = rep(NA, max_iter)
  rho03value = rep(NA, max_iter)
  rho04value = rep(NA, max_iter)
  rho05value = rep(NA, max_iter)
  rho06value = rep(NA, max_iter)
  rho07value = rep(NA, max_iter)
  
  rho1value = rep(NA, max_iter)
  rho2value = rep(NA, max_iter)
  rho3value = rep(NA, max_iter)
  rho4value = rep(NA, max_iter)
  rho5value = rep(NA, max_iter)
  rho6value = rep(NA, max_iter)
  rho7value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
  r_primal04value = rep(NA, max_iter) ; s_dual04value = rep(NA, max_iter)
  r_primal05value = rep(NA, max_iter) ; s_dual05value = rep(NA, max_iter)
  r_primal06value = rep(NA, max_iter) ; s_dual06value = rep(NA, max_iter)
  r_primal07value = rep(NA, max_iter) ; s_dual07value = rep(NA, max_iter)
  
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
  r_primal4value = rep(NA, max_iter) ; s_dual4value = rep(NA, max_iter)
  r_primal5value = rep(NA, max_iter) ; s_dual5value = rep(NA, max_iter)
  r_primal6value = rep(NA, max_iter) ; s_dual6value = rep(NA, max_iter)
  r_primal7value = rep(NA, max_iter) ; s_dual7value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  B2 <- matrix(0, p, 1)
  B3 <- matrix(0, p, 1)
  B4 <- matrix(0, p, 1)
  B5 <- matrix(0, p, 1)
  B6 <- matrix(0, p, 1)
  B7 <- matrix(0, p, 1)
  
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma02 <- matrix(0, p, 1)
  Gamma03 <- matrix(0, p, 1)
  Gamma04 <- matrix(0, p, 1)
  Gamma05 <- matrix(0, p, 1)
  Gamma06 <- matrix(0, p, 1)
  Gamma07 <- matrix(0, p, 1)
  
  Gamma1 <- matrix(0, p, 1)
  Gamma2 <- matrix(0, p, 1)
  Gamma3 <- matrix(0, p, 1)
  Gamma4 <- matrix(0, p, 1)
  Gamma5 <- matrix(0, p, 1)
  Gamma6 <- matrix(0, p, 1)
  Gamma7 <- matrix(0, p, 1)
  
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  U2 <- matrix(0, p, 1)
  U3 <- matrix(0, p, 1)
  U4 <- matrix(0, p, 1)
  U5 <- matrix(0, p, 1)
  U6 <- matrix(0, p, 1)
  U7 <- matrix(0, p, 1)
  
  V1 <- matrix(0, p, 1)
  V2 <- matrix(0, p, 1)
  V3 <- matrix(0, p, 1)
  V4 <- matrix(0, p, 1)
  V5 <- matrix(0, p, 1)
  V6 <- matrix(0, p, 1)
  V7 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0) ; X1tX1 = crossprod(X1,X1) ; X2tX2 = crossprod(X2,X2) ; X3tX3 = crossprod(X3,X3) ; X4tX4 = crossprod(X4,X4) ; X5tX5 = crossprod(X5,X5) ; X6tX6 = crossprod(X6,X6) ; X7tX7 = crossprod(X7,X7)
    X0ty0 = crossprod(X0,y0) ; X1ty1 = crossprod(X1,y1) ; X2ty2 = crossprod(X2,y2) ; X3ty3 = crossprod(X3,y3) ; X4ty4 = crossprod(X4,y4) ; X5ty5 = crossprod(X5,y5) ; X6ty6 = crossprod(X6,y6) ; X7ty7 = crossprod(X7,y7)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    # if(iter %% 100 == 0){cat('iter = ', iter, '\n') }
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma03_old = Gamma03
    Gamma04_old = Gamma04
    Gamma05_old = Gamma05
    Gamma06_old = Gamma06
    Gamma07_old = Gamma07
    
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    Gamma3_old = Gamma3
    Gamma4_old = Gamma4
    Gamma5_old = Gamma5
    Gamma6_old = Gamma6
    Gamma7_old = Gamma7
    
    #### Step 1: B0 update ####
    lhs = 1/N*X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06 + rho07)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6) + rho07 * (Gamma07 - U7)
    B0 <- solve(lhs, rhs)
    
    #### Step 2: B1 update ####
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### Step 3: B2 update ####
    lhs = 1/N*X2tX2 + rho2*Ip
    rhs = 1/N*X2ty2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### Step 4: B3 update ####
    lhs = 1/N*X3tX3 + rho3*Ip
    rhs = 1/N*X3ty3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### Step 5: B4 update ####
    lhs = 1/N*X4tX4 + rho4*Ip
    rhs = 1/N*X4ty4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### Step 6: B5 update ####
    lhs = 1/N*X5tX5 + rho5*Ip
    rhs = 1/N*X5ty5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    ####  B6 update ####
    lhs = 1/N*X6tX6 + rho6*Ip
    rhs = 1/N*X6ty6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    
    ####  B7 update ####
    lhs = 1/N*X7tX7 + rho7*Ip
    rhs = 1/N*X7ty7 + rho7*(Gamma7 - V7)
    B7 <- solve(lhs, rhs)
    
    
    #### Step 7: Gamma00 update ####
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Step 8: Gamma01 update ####
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    
    #### Step 9: Gamma02 update ####
    Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
    
    #### Step 10: Gamma03 update ####
    Gamma03 = Gamma3 + soft_threshold(B0+U3-Gamma3, lambda3/rho03)
    
    #### Step 11: Gamma04 update ####
    Gamma04 = Gamma4 + soft_threshold(B0+U4-Gamma4, lambda4/rho04)
    
    #### Step 12: Gamma05 update ####
    Gamma05 = Gamma5 + soft_threshold(B0+U5-Gamma5, lambda5/rho05)
    
    #### Gamma06 update ####
    Gamma06 = Gamma6 + soft_threshold(B0+U6-Gamma6, lambda6/rho06)
    
    #### Gamma07 update ####
    Gamma07 = Gamma7 + soft_threshold(B0+U7-Gamma7, lambda7/rho07)
    
    
    
    #### Step 13: Gamma1 update ####
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
    
    #### Step 14: Gamma2 update ####
    Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
    
    #### Step 15: Gamma3 update ####
    Gamma3 = Gamma03 + soft_threshold(B3+V3-Gamma03, lambda3/rho3)
    
    #### Step 16: Gamma4 update ####
    Gamma4 = Gamma04 + soft_threshold(B4+V4-Gamma04, lambda4/rho4)
    
    #### Step 17: Gamma5 update ####
    Gamma5 = Gamma05 + soft_threshold(B5+V5-Gamma05, lambda5/rho5)
    
    #### Step  Gamma6 update ####
    Gamma6 = Gamma06 + soft_threshold(B6+V6-Gamma06, lambda6/rho6)
    
    #### Step  Gamma7 update ####
    Gamma7 = Gamma07 + soft_threshold(B7+V7-Gamma07, lambda7/rho7)
    
    #### Step 18–28: Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    U3 <- U3 + B0 - Gamma03
    U4 <- U4 + B0 - Gamma04
    U5 <- U5 + B0 - Gamma05
    U6 <- U6 + B0 - Gamma06
    U7 <- U7 + B0 - Gamma07
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    V3 <- V3 + B3 - Gamma3
    V4 <- V4 + B4 - Gamma4
    V5 <- V5 + B5 - Gamma5
    V6 <- V6 + B6 - Gamma6
    V7 <- V7 + B7 - Gamma7
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
    r_primal03 = l2norm(B0 - Gamma03) ; s_dual03 = rho03 * l2norm(Gamma03 - Gamma03_old)
    r_primal04 = l2norm(B0 - Gamma04) ; s_dual04 = rho04 * l2norm(Gamma04 - Gamma04_old)
    r_primal05 = l2norm(B0 - Gamma05) ; s_dual05 = rho05 * l2norm(Gamma05 - Gamma05_old)
    r_primal06 = l2norm(B0 - Gamma06) ; s_dual06 = rho06 * l2norm(Gamma06 - Gamma06_old)
    r_primal07 = l2norm(B0 - Gamma07) ; s_dual07 = rho07 * l2norm(Gamma07 - Gamma07_old)
    
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
    r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
    r_primal3 = l2norm(B3 - Gamma3) ; s_dual3 = rho3 * l2norm(Gamma3 - Gamma3_old)
    r_primal4 = l2norm(B4 - Gamma4) ; s_dual4 = rho4 * l2norm(Gamma4 - Gamma4_old)
    r_primal5 = l2norm(B5 - Gamma5) ; s_dual5 = rho5 * l2norm(Gamma5 - Gamma5_old)
    r_primal6 = l2norm(B6 - Gamma6) ; s_dual6 = rho6 * l2norm(Gamma6 - Gamma6_old)
    r_primal7 = l2norm(B7 - Gamma7) ; s_dual7 = rho7 * l2norm(Gamma7 - Gamma7_old)
    
    # if (verbose && iter %% 1000 == 0) {
    #   cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid02: %.2f | Dual Resid02: %.2f | Primal Resid03: %.2f |Dual Resid03: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f | Primal Resid2: %.2f | Dual Resid2: %.2f | Primal Resid3: %.2f | Dual Resid3: %.2f \n",
    #               iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal03,base=10), log(s_dual03,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10), log(r_primal3,base=10), log(s_dual3,base=10)))
    # }
    
    # objvalue[iter] = 1/(2*N)*(l2norm_sq(y0 - X0 %*% B0) + l2norm_sq(y1 - X1 %*% B1) + l2norm_sq(y2 - X2 %*% B2) + l2norm_sq(y3 - X3 %*% B3) + l2norm_sq(y4 - X4 %*% B4) + l2norm_sq(y5 - X5 %*% B5)) +
    #   lambda0 * sum(abs(Gamma00)) + lambda1 * sum(abs(Gamma1-Gamma00)) + lambda2 * sum(abs(Gamma2-Gamma00)) + lambda3 * sum(abs(Gamma3-Gamma00)) + lambda4 * sum(abs(Gamma4-Gamma00)) + lambda5 * sum(abs(Gamma5-Gamma00))  
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
    r_primal03value[iter] = r_primal03 ; s_dual03value[iter] = s_dual03
    r_primal04value[iter] = r_primal04 ; s_dual04value[iter] = s_dual04
    r_primal05value[iter] = r_primal05 ; s_dual05value[iter] = s_dual05
    r_primal06value[iter] = r_primal06 ; s_dual06value[iter] = s_dual06
    r_primal07value[iter] = r_primal07 ; s_dual07value[iter] = s_dual07
    
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
    r_primal3value[iter] = r_primal3 ; s_dual3value[iter] = s_dual3
    r_primal4value[iter] = r_primal4 ; s_dual4value[iter] = s_dual4
    r_primal5value[iter] = r_primal5 ; s_dual5value[iter] = s_dual5
    r_primal6value[iter] = r_primal6 ; s_dual6value[iter] = s_dual6
    r_primal7value[iter] = r_primal7 ; s_dual7value[iter] = s_dual7
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal02 > mu * s_dual02){
      rho02 = min(rho02 * tau_incr, rho_max)
      U2 = U2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual02 > mu * r_primal02){
      rho02 = rho02 / tau_decr
      U2 = U2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal03 > mu * s_dual03){
      rho03 = min(rho03 * tau_incr, rho_max)
      U3 = U3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual03 > mu * r_primal03){
      rho03 = rho03 / tau_decr
      U3 = U3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal04 > mu * s_dual04){
      rho04 = min(rho04 * tau_incr, rho_max)
      U4 = U4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual04 > mu * r_primal04){
      rho04 = rho04 / tau_decr
      U4 = U4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal05 > mu * s_dual05){
      rho05 = min(rho05 * tau_incr, rho_max)
      U5 = U5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual05 > mu * r_primal05){
      rho05 = rho05 / tau_decr
      U5 = U5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal06 > mu * s_dual06){
      rho06 = min(rho06 * tau_incr, rho_max)
      U6 = U6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual06 > mu * r_primal06){
      rho06 = rho06 / tau_decr
      U6 = U6 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal07 > mu * s_dual07){
      rho07 = min(rho07 * tau_incr, rho_max)
      U7 = U7 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual07 > mu * r_primal07){
      rho07 = rho07 / tau_decr
      U7 = U7 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal2 > mu * s_dual2){
      rho2 = min(rho2 * tau_incr, rho_max)
      V2 = V2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual2 > mu * r_primal2){
      rho2 = rho2 / tau_decr
      V2 = V2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal3 > mu * s_dual3){
      rho3 = min(rho3 * tau_incr, rho_max)
      V3 = V3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual3 > mu * r_primal3){
      rho3 = rho3 / tau_decr
      V3 = V3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal4 > mu * s_dual4){
      rho4 = min(rho4 * tau_incr, rho_max)
      V4 = V4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual4 > mu * r_primal4){
      rho4 = rho4 / tau_decr
      V4 = V4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal5 > mu * s_dual5){
      rho5 = min(rho5 * tau_incr, rho_max)
      V5 = V5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual5 > mu * r_primal5){
      rho5 = rho5 / tau_decr
      V5 = V5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal6 > mu * s_dual6){
      rho6 = min(rho6 * tau_incr, rho_max)
      V6 = V6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual6 > mu * r_primal6){
      rho6 = rho6 / tau_decr
      V6 = V6 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal7 > mu * s_dual7){
      rho7 = min(rho7 * tau_incr, rho_max)
      V7 = V7 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual7 > mu * r_primal7){
      rho7 = rho7 / tau_decr
      V7 = V7 * tau_decr # change of variable required for the scaled variable
    }
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho02value[iter] = rho02
    rho03value[iter] = rho03
    rho04value[iter] = rho04
    rho05value[iter] = rho05
    rho06value[iter] = rho06
    rho07value[iter] = rho07
    
    rho1value[iter] = rho1
    rho2value[iter] = rho2
    rho3value[iter] = rho3
    rho4value[iter] = rho4
    rho5value[iter] = rho5
    rho6value[iter] = rho6
    rho7value[iter] = rho7
    
    if (max(r_primal00, 
            r_primal01, r_primal02, r_primal03, r_primal04, r_primal05, r_primal06, r_primal07, 
            r_primal1, r_primal2, r_primal3, r_primal4, r_primal5, r_primal6, r_primal7) < tol_prim 
        &&
        max(s_dual00, 
            s_dual01, s_dual02, s_dual03, s_dual04, s_dual05, s_dual06, s_dual07, 
            s_dual1, s_dual2, s_dual3, s_dual4, s_dual5, s_dual6, s_dual7) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6, B7 = B7,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06, Gamma07 = Gamma07,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6, Gamma7 = Gamma7,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6, U7 = U7,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6, V7 = V7,
              #objvalue = objvalue, 
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, rho04 = rho04value, rho05 = rho05value, rho06 = rho06value, rho07 = rho07value,
              rho1 = rho1value, rho2 = rho2value, rho3 = rho3value, rho4 = rho4value, rho5 = rho5value, rho6 = rho6value, rho7 = rho7value,
              
              r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              r_primal04 = r_primal04value, s_dual04 = s_dual04value,
              r_primal05 = r_primal05value, s_dual05 = s_dual05value,
              r_primal06 = r_primal06value, s_dual06 = s_dual06value,
              r_primal07 = r_primal07value, s_dual07 = s_dual07value,
              
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value,
              r_primal4 = r_primal4value, s_dual4 = s_dual4value,
              r_primal5 = r_primal5value, s_dual5 = s_dual5value,
              r_primal6 = r_primal6value, s_dual6 = s_dual6value,
              r_primal7 = r_primal7value, s_dual7 = s_dual7value))
}

# cvx_stl_diffweight = function(X0,y0,X1,y1,
#                               lambda0, lambda1){ # CVXR implementation of He et al. 2024 method with different weight for source & targe
#   
#   n_target = nrow(X0) ; n_source = nrow(X1)
#   p = ncol(X0)
#   
#   w_source = Variable(p,1) ; w_target = Variable(p,1)
#   loss = 1/(2*n_target)*(sum_squares(X0 %*% w_target - y0)) + 1/(2*n_source)*(sum_squares(X1 %*% w_source - y1))
#   objective = Minimize(loss + lambda0 * p_norm(w_target, p = 1) + lambda1 * p_norm(w_target - w_source, p = 1))
#   problem = Problem(objective)
#   result = solve(problem)
#   w_opt_target = as.matrix(result$getValue(w_target))
#   w_opt_source = as.matrix(result$getValue(w_source))
#   
#   return(list(B0 = w_opt_target, B1 = w_opt_source, result = result))
# }
# 
# admm_stl_diffweight <- function(X0, y0, X1, y1, 
#                                 X0tX0 = NULL, X0ty0 = NULL, X1tX1 = NULL, X1ty1 = NULL,
#                                 lambda0, lambda1, 
#                                 rho00 = 1, rho01 = 1, rho1 = 1, 
#                                 mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
#                                 max_iter = 2000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){
#   
#   l2norm = function(x){sqrt(sum(x^2))}
#   l2norm_sq = function(x){sum(x^2)}
#   
#   n_target = nrow(X0)
#   n_source = nrow(X1)
#   
#   objvalue = rep(NA, max_iter)
#   rho00value = rep(NA, max_iter)
#   rho01value = rep(NA, max_iter)
#   rho1value = rep(NA, max_iter)
#   
#   r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
#   r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
#   r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
#   
#   p <- ncol(X0) # same number of features across feature and target
#   
#   B0 <- matrix(0, p, 1)
#   B1 <- matrix(0, p, 1)
#   Gamma00 <- matrix(0, p, 1)
#   Gamma01 <- matrix(0, p, 1)
#   Gamma1 <- matrix(0, p, 1)
#   U0 <- matrix(0, p, 1)
#   U1 <- matrix(0, p, 1)
#   V1 <- matrix(0, p, 1)
#   
#   if(is.null(X0tX0)){
#     X0tX0 = crossprod(X0,X0) ; X1tX1 = crossprod(X1,X1)
#     X0ty0 = crossprod(X0,y0) ; X1ty1 = crossprod(X1,y1)
#   }
#   Ip = diag(p)
#   
#   for (iter in 1:max_iter) {
#     
#     Gamma00_old = Gamma00
#     Gamma01_old = Gamma01
#     Gamma1_old = Gamma1
#     
#     #### Step 1: B0 update ####
#     lhs = 1/n_target*X0tX0 + (rho00 + rho01)*Ip
#     rhs = 1/n_target*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1)
#     B0 <- solve(lhs, rhs)
#     
#     #### Step 2: B1 update ####
#     lhs = 1/n_source*X1tX1 + rho1*Ip
#     rhs = 1/n_source*X1ty1 + rho1*(Gamma1 - V1)
#     B1 <- solve(lhs, rhs)
#     
#     #### Step 3: Gamma00 update ####
#     Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
#     
#     #### Step 4: Gamma01 update ####
#     Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
#     
#     #### Step 5: Gamma1 update ####
#     Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
#     
#     #### Step 6–8: Dual updates ####
#     U0 <- U0 + B0 - Gamma00
#     U1 <- U1 + B0 - Gamma01
#     V1 <- V1 + B1 - Gamma1
#     
#     #### Convergence check: primal & dual ####
#     r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
#     r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
#     r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
#     
#     if (verbose && iter %% 500 == 0) {
#       cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f \n",
#                   iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal1,base=10), log(s_dual1,base=10) ))
#     }
#     
#     
#     objvalue[iter] = 1/(2*n_target)*l2norm_sq(y0 - X0 %*% Gamma00) + 1/(2*n_source)*l2norm_sq(y1 - X1 %*% Gamma1) + lambda0 * sum(abs(Gamma00)) + lambda1 * sum(abs(Gamma1-Gamma00)) # objective value
#     
#     r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
#     r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
#     r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
#     
#     #Adaptive rho adjustment (Boyd 2011)
#     if(r_primal00 > mu * s_dual00){
#       rho00 = min(rho00 * tau_incr, rho_max)
#       U0 = U0 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual00 > mu * r_primal00){
#       rho00 = rho00 / tau_decr
#       U0 = U0 * tau_decr # change of variable required for the scaled variable
#     }
#     
#     if(r_primal01 > mu * s_dual01){
#       rho01 = min(rho01 * tau_incr, rho_max)
#       U1 = U1 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual01 > mu * r_primal01){
#       rho01 = rho01 / tau_decr
#       U1 = U1 * tau_decr # change of variable required for the scaled variable
#     }
#     
#     if(r_primal1 > mu * s_dual1){
#       rho1 = min(rho1 * tau_incr, rho_max)
#       V1 = V1 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual1 > mu * r_primal1){
#       rho1 = rho1 / tau_decr
#       V1 = V1 * tau_decr # change of variable required for the scaled variable
#     }
#     
#     
#     rho00value[iter] = rho00
#     rho01value[iter] = rho01
#     rho1value[iter] = rho1
#     
#     if (max(r_primal00, r_primal01, r_primal1) < tol_prim && max(s_dual00, s_dual01, s_dual1) < tol_dual) {
#       if(verbose){
#         cat("ADMM converged at iteration", iter, "\n")
#       }
#       break
#     }
#   }
#   
#   return(list(B0 = B0, B1 = B1, 
#               Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma1 = Gamma1,
#               U0 = U0, U1 = U1, V1 = V1, 
#               objvalue = objvalue, 
#               rho00 = rho00value, rho01 = rho01value, rho1 = rho1value,
#               r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
#               r_primal01 = r_primal01value, s_dual01 = s_dual01value,
#               r_primal1 = r_primal1value, s_dual1 = s_dual1value))
# }
# 
# admm_stl_twosource_diffweight <- function(X0, y0, X1, y1, X2, y2,
#                                           X0tX0 = NULL, X0ty0 = NULL, X1tX1 = NULL, X1ty1 = NULL, X2tX2 = NULL, X2ty2 = NULL,
#                                           lambda0, lambda1, lambda2,
#                                           rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1,
#                                           mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
#                                           max_iter = 2000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){
#   
#   l2norm = function(x){sqrt(sum(x^2))}
#   l2norm_sq = function(x){sum(x^2)}
#   
#   n_target = nrow(X0)
#   n_source = nrow(X1)
#   
#   objvalue = rep(NA, max_iter)
#   
#   rho00value = rep(NA, max_iter)
#   rho01value = rep(NA, max_iter)
#   rho02value = rep(NA, max_iter)
#   
#   rho1value = rep(NA, max_iter)
#   rho2value = rep(NA, max_iter)
#   
#   r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
#   r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
#   r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
#   
#   r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
#   r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
#   
#   p <- ncol(X0)
#   
#   B0 <- matrix(0, p, 1)
#   B1 <- matrix(0, p, 1)
#   B2 <- matrix(0, p, 1)
#   
#   Gamma00 <- matrix(0, p, 1)
#   Gamma01 <- matrix(0, p, 1)
#   Gamma02 <- matrix(0, p, 1)
#   
#   Gamma1 <- matrix(0, p, 1)
#   Gamma2 <- matrix(0, p, 1)
#   
#   U0 <- matrix(0, p, 1)
#   U1 <- matrix(0, p, 1)
#   U2 <- matrix(0, p, 1)
#   
#   V1 <- matrix(0, p, 1)
#   V2 <- matrix(0, p, 1)
#   
#   if(is.null(X0tX0)){
#     X0tX0 = crossprod(X0,X0) ; X1tX1 = crossprod(X1,X1) ; X2tX2 = crossprod(X2,X2) 
#     X0ty0 = crossprod(X0,y0) ; X1ty1 = crossprod(X1,y1) ; X2ty2 = crossprod(X2,y2) 
#   }
#   Ip = diag(p)
#   
#   for (iter in 1:max_iter) {
#     
#     Gamma00_old = Gamma00
#     Gamma01_old = Gamma01
#     Gamma02_old = Gamma02
#     
#     Gamma1_old = Gamma1
#     Gamma2_old = Gamma2
#     
#     #### Step 1: B0 update ####
#     lhs = 1/n_target*X0tX0 + (rho00 + rho01 + rho02)*Ip
#     rhs = 1/n_target*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2)
#     B0 <- solve(lhs, rhs)
#     
#     #### Step 2: B1 update ####
#     lhs = 1/n_source*X1tX1 + rho1*Ip
#     rhs = 1/n_source*X1ty1 + rho1*(Gamma1 - V1)
#     B1 <- solve(lhs, rhs)
#     
#     #### Step 3: B2 update ####
#     lhs = 1/n_source*X2tX2 + rho2*Ip
#     rhs = 1/n_source*X2ty2 + rho2*(Gamma2 - V2)
#     B2 <- solve(lhs, rhs)
#     
#     #### Step 4: Gamma00 update ####
#     Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
#     
#     #### Step 5: Gamma01 update ####
#     Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
#     
#     #### Step 6: Gamma02 update ####
#     Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
#     
#     #### Step 7: Gamma1 update ####
#     Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
#     
#     #### Step 8: Gamma2 update ####
#     Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
#     
#     #### Step 9–13: Dual updates ####
#     U0 <- U0 + B0 - Gamma00
#     U1 <- U1 + B0 - Gamma01
#     U2 <- U2 + B0 - Gamma02
#     
#     V1 <- V1 + B1 - Gamma1
#     V2 <- V2 + B2 - Gamma2
#     
#     #### Convergence check: primal & dual ####
#     r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
#     r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
#     r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
#     
#     r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
#     r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
#     
#     if (verbose && iter %% 500 == 0) {
#       cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid02: %.2f | Dual Resid02: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f | Primal Resid2: %.2f | Dual Resid2: %.2f \n",
#                   iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10)))
#     }
#     objvalue[iter] = 1/(2*n_target)*l2norm_sq(y0 - X0 %*% Gamma00) + 1/(2*n_source)*l2norm_sq(y1 - X1 %*% Gamma1) + 1/(2*n_source)*l2norm_sq(y2 - X2 %*% Gamma2) + lambda0 * sum(abs(Gamma00)) + lambda1 * sum(abs(Gamma1-Gamma00)) + lambda2 * sum(abs(Gamma2-Gamma00)) # objective value
#     
#     r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
#     r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
#     r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
#     
#     r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
#     r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
#     
#     #Adaptive rho adjustment (Boyd 2011)
#     if(r_primal00 > mu * s_dual00){
#       rho00 = min(rho00 * tau_incr, rho_max)
#       U0 = U0 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual00 > mu * r_primal00){
#       rho00 = rho00 / tau_decr
#       U0 = U0 * tau_decr # change of variable required for the scaled variable
#     }
#     
#     if(r_primal01 > mu * s_dual01){
#       rho01 = min(rho01 * tau_incr, rho_max)
#       U1 = U1 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual01 > mu * r_primal01){
#       rho01 = rho01 / tau_decr
#       U1 = U1 * tau_decr # change of variable required for the scaled variable
#     }
#     if(r_primal02 > mu * s_dual02){
#       rho02 = min(rho02 * tau_incr, rho_max)
#       U2 = U2 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual02 > mu * r_primal02){
#       rho02 = rho02 / tau_decr
#       U2 = U2 * tau_decr # change of variable required for the scaled variable
#     }
#     
#     
#     if(r_primal1 > mu * s_dual1){
#       rho1 = min(rho1 * tau_incr, rho_max)
#       V1 = V1 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual1 > mu * r_primal1){
#       rho1 = rho1 / tau_decr
#       V1 = V1 * tau_decr # change of variable required for the scaled variable
#     }
#     if(r_primal2 > mu * s_dual2){
#       rho2 = min(rho2 * tau_incr, rho_max)
#       V2 = V2 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual2 > mu * r_primal2){
#       rho2 = rho2 / tau_decr
#       V2 = V2 * tau_decr # change of variable required for the scaled variable
#     }
#     
#     
#     rho00value[iter] = rho00
#     rho01value[iter] = rho01
#     rho02value[iter] = rho02
#     
#     rho1value[iter] = rho1
#     rho2value[iter] = rho2
#     
#     if (max(r_primal00, r_primal01, r_primal02, r_primal1, r_primal2) < tol_prim && max(s_dual00, s_dual01, s_dual02, s_dual1, s_dual2) < tol_dual) {
#       if(verbose){
#         cat("ADMM converged at iteration", iter, "\n")
#       }
#       break
#     }
#   }
#   
#   return(list(B0 = B0, B1 = B1, B2 = B2,
#               Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02,
#               Gamma1 = Gamma1, Gamma2 = Gamma2,
#               U0 = U0, U1 = U1, U2 = U2,
#               V1 = V1, V2 = V2,
#               objvalue = objvalue, 
#               rho00 = rho00value, rho01 = rho01value, rho02 = rho02value,
#               rho1 = rho1value, rho2 = rho2value,
#               
#               r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
#               r_primal01 = r_primal01value, s_dual01 = s_dual01value,
#               r_primal02 = r_primal02value, s_dual02 = s_dual02value,
#               
#               r_primal1 = r_primal1value, s_dual1 = s_dual1value,
#               r_primal2 = r_primal2value, s_dual2 = s_dual2value))
# }
# 
# admm_stl_threesource_diffweight <- function(X0, y0, X1, y1, X2, y2, X3, y3,
#                                             X0tX0 = NULL, X0ty0 = NULL, X1tX1 = NULL, X1ty1 = NULL, X2tX2 = NULL, X2ty2 = NULL, X3tX3 = NULL, X3ty3 = NULL,
#                                             lambda0, lambda1, lambda2, lambda3,
#                                             rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho1 = 1, rho2 = 1, rho3 = 1,
#                                             mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
#                                             max_iter = 5000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){ 
#   
#   l2norm = function(x){sqrt(sum(x^2))}
#   l2norm_sq = function(x){sum(x^2)}
#   
#   n_target = nrow(X0)
#   n_source = nrow(X1)
#   
#   objvalue = rep(NA, max_iter)
#   
#   rho00value = rep(NA, max_iter)
#   rho01value = rep(NA, max_iter)
#   rho02value = rep(NA, max_iter)
#   rho03value = rep(NA, max_iter)
#   
#   rho1value = rep(NA, max_iter)
#   rho2value = rep(NA, max_iter)
#   rho3value = rep(NA, max_iter)
#   
#   r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
#   r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
#   r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
#   r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
#   
#   r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
#   r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
#   r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
#   
#   p <- ncol(X0)
#   
#   B0 <- matrix(0, p, 1)
#   B1 <- matrix(0, p, 1)
#   B2 <- matrix(0, p, 1)
#   B3 <- matrix(0, p, 1)
#   
#   Gamma00 <- matrix(0, p, 1)
#   Gamma01 <- matrix(0, p, 1)
#   Gamma02 <- matrix(0, p, 1)
#   Gamma03 <- matrix(0, p, 1)
#   
#   Gamma1 <- matrix(0, p, 1)
#   Gamma2 <- matrix(0, p, 1)
#   Gamma3 <- matrix(0, p, 1)
#   
#   U0 <- matrix(0, p, 1)
#   U1 <- matrix(0, p, 1)
#   U2 <- matrix(0, p, 1)
#   U3 <- matrix(0, p, 1)
#   
#   V1 <- matrix(0, p, 1)
#   V2 <- matrix(0, p, 1)
#   V3 <- matrix(0, p, 1)
#   
#   if(is.null(X0tX0)){
#     X0tX0 = crossprod(X0,X0) ; X1tX1 = crossprod(X1,X1) ; X2tX2 = crossprod(X2,X2) ; X3tX3 = crossprod(X3,X3)
#     X0ty0 = crossprod(X0,y0) ; X1ty1 = crossprod(X1,y1) ; X2ty2 = crossprod(X2,y2) ; X3ty3 = crossprod(X3,y3)
#   }
#   Ip = diag(p)
#   
#   for (iter in 1:max_iter) {
#     
#     Gamma00_old = Gamma00
#     Gamma01_old = Gamma01
#     Gamma02_old = Gamma02
#     Gamma03_old = Gamma03
#     
#     Gamma1_old = Gamma1
#     Gamma2_old = Gamma2
#     Gamma3_old = Gamma3
#     
#     #### Step 1: B0 update ####
#     lhs = 1/n_target*X0tX0 + (rho00 + rho01 + rho02 + rho03)*Ip
#     rhs = 1/n_target*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3)
#     B0 <- solve(lhs, rhs)
#     
#     #### Step 2: B1 update ####
#     lhs = 1/n_source*X1tX1 + rho1*Ip
#     rhs = 1/n_source*X1ty1 + rho1*(Gamma1 - V1)
#     B1 <- solve(lhs, rhs)
#     
#     #### Step 3: B2 update ####
#     lhs = 1/n_source*X2tX2 + rho2*Ip
#     rhs = 1/n_source*X2ty2 + rho2*(Gamma2 - V2)
#     B2 <- solve(lhs, rhs)
#     
#     #### Step 4: B3 update ####
#     lhs = 1/n_source*X3tX3 + rho3*Ip
#     rhs = 1/n_source*X3ty3 + rho3*(Gamma3 - V3)
#     B3 <- solve(lhs, rhs)
#     
#     #### Step 5: Gamma00 update ####
#     Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
#     
#     #### Step 6: Gamma01 update ####
#     Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
#     
#     #### Step 7: Gamma02 update ####
#     Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
#     
#     #### Step 8: Gamma03 update ####
#     Gamma03 = Gamma3 + soft_threshold(B0+U3-Gamma3, lambda3/rho03)
#     
#     #### Step 9: Gamma1 update ####
#     Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
#     
#     #### Step 10: Gamma2 update ####
#     Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
#     
#     #### Step 11: Gamma3 update ####
#     Gamma3 = Gamma03 + soft_threshold(B3+V3-Gamma03, lambda3/rho3)
#     
#     #### Step 12–18: Dual updates ####
#     U0 <- U0 + B0 - Gamma00
#     U1 <- U1 + B0 - Gamma01
#     U2 <- U2 + B0 - Gamma02
#     U3 <- U3 + B0 - Gamma03
#     
#     V1 <- V1 + B1 - Gamma1
#     V2 <- V2 + B2 - Gamma2
#     V3 <- V3 + B3 - Gamma3
#     
#     #### Convergence check: primal & dual ####
#     r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
#     r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
#     r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
#     r_primal03 = l2norm(B0 - Gamma03) ; s_dual03 = rho03 * l2norm(Gamma03 - Gamma03_old)
#     
#     r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
#     r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
#     r_primal3 = l2norm(B3 - Gamma3) ; s_dual3 = rho3 * l2norm(Gamma3 - Gamma3_old)
#     
#     if (verbose && iter %% 500 == 0) {
#       cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid02: %.2f | Dual Resid02: %.2f | Primal Resid03: %.2f |Dual Resid03: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f | Primal Resid2: %.2f | Dual Resid2: %.2f | Primal Resid3: %.2f | Dual Resid3: %.2f \n",
#                   iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal03,base=10), log(s_dual03,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10), log(r_primal3,base=10), log(s_dual3,base=10)))
#     }
#     
#     objvalue[iter] = 1/(2*n_target)*l2norm_sq(y0 - X0 %*% Gamma00) + 1/(2*n_source)*l2norm_sq(y1 - X1 %*% Gamma1) + 1/(2*n_source)*l2norm_sq(y2 - X2 %*% Gamma2) + 1/(2*n_source)*l2norm_sq(y3 - X3 %*% Gamma3) + lambda0 * sum(abs(Gamma00)) + lambda1 * sum(abs(Gamma1-Gamma00)) + lambda2 * sum(abs(Gamma2-Gamma00)) + lambda3 * sum(abs(Gamma3-Gamma00)) # objective value
#     
#     r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
#     r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
#     r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
#     r_primal03value[iter] = r_primal03 ; s_dual03value[iter] = s_dual03
#     
#     r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
#     r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
#     r_primal3value[iter] = r_primal3 ; s_dual3value[iter] = s_dual3
#     
#     #Adaptive rho adjustment (Boyd 2011)
#     if(r_primal00 > mu * s_dual00){
#       rho00 = min(rho00 * tau_incr, rho_max)
#       U0 = U0 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual00 > mu * r_primal00){
#       rho00 = rho00 / tau_decr
#       U0 = U0 * tau_decr # change of variable required for the scaled variable
#     }
#     
#     if(r_primal01 > mu * s_dual01){
#       rho01 = min(rho01 * tau_incr, rho_max)
#       U1 = U1 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual01 > mu * r_primal01){
#       rho01 = rho01 / tau_decr
#       U1 = U1 * tau_decr # change of variable required for the scaled variable
#     }
#     if(r_primal02 > mu * s_dual02){
#       rho02 = min(rho02 * tau_incr, rho_max)
#       U2 = U2 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual02 > mu * r_primal02){
#       rho02 = rho02 / tau_decr
#       U2 = U2 * tau_decr # change of variable required for the scaled variable
#     }
#     if(r_primal03 > mu * s_dual03){
#       rho03 = min(rho03 * tau_incr, rho_max)
#       U3 = U3 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual03 > mu * r_primal03){
#       rho03 = rho03 / tau_decr
#       U3 = U3 * tau_decr # change of variable required for the scaled variable
#     }
#     
#     if(r_primal1 > mu * s_dual1){
#       rho1 = min(rho1 * tau_incr, rho_max)
#       V1 = V1 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual1 > mu * r_primal1){
#       rho1 = rho1 / tau_decr
#       V1 = V1 * tau_decr # change of variable required for the scaled variable
#     }
#     if(r_primal2 > mu * s_dual2){
#       rho2 = min(rho2 * tau_incr, rho_max)
#       V2 = V2 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual2 > mu * r_primal2){
#       rho2 = rho2 / tau_decr
#       V2 = V2 * tau_decr # change of variable required for the scaled variable
#     }
#     if(r_primal3 > mu * s_dual3){
#       rho3 = min(rho3 * tau_incr, rho_max)
#       V3 = V3 / tau_incr # change of variable required for the scaled variable
#     }else if(s_dual3 > mu * r_primal3){
#       rho3 = rho3 / tau_decr
#       V3 = V3 * tau_decr # change of variable required for the scaled variable
#     }
#     
#     
#     rho00value[iter] = rho00
#     rho01value[iter] = rho01
#     rho02value[iter] = rho02
#     rho03value[iter] = rho03
#     
#     rho1value[iter] = rho1
#     rho2value[iter] = rho2
#     rho3value[iter] = rho3
#     
#     if (max(r_primal00, r_primal01, r_primal02, r_primal03, r_primal1, r_primal2, r_primal3) < tol_prim && max(s_dual00, s_dual01, s_dual02, s_dual03, s_dual1, s_dual2, s_dual3) < tol_dual) {
#       if(verbose){
#         cat("ADMM converged at iteration", iter, "\n")
#       }
#       break
#     }
#   }
#   
#   return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, 
#               Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, 
#               Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3,
#               U0 = U0, U1 = U1, U2 = U2, U3 = U3, 
#               V1 = V1, V2 = V2, V3 = V3, 
#               objvalue = objvalue, 
#               rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, 
#               rho1 = rho1value, rho2 = rho2value, rho3 = rho3value,
#               
#               r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
#               r_primal01 = r_primal01value, s_dual01 = s_dual01value,
#               r_primal02 = r_primal02value, s_dual02 = s_dual02value,
#               r_primal03 = r_primal03value, s_dual03 = s_dual03value,
#               
#               r_primal1 = r_primal1value, s_dual1 = s_dual1value,
#               r_primal2 = r_primal2value, s_dual2 = s_dual2value,
#               r_primal3 = r_primal3value, s_dual3 = s_dual3value))
# }






# #### Diagnosis ####
# 
# load("~/Library/CloudStorage/GoogleDrive-96limtotoro@gmail.com/My Drive/research_multitask/STL_check/p50/ntarget100_nsource100/onesource.Rdata")
# 
# # data generation 
# set.seed(1)
# sigma = 0.2
# lambda0 = 0.1 ; lambda1 = 0.1
# 
# X_target = matrix(rnorm(n_target * p, mean = 0, sd = 1), n_target, p) ; X_target = scale(X_target, center = T, scale = F) #; X_target_center = attr(X_target, "scaled:center")
# Y_target = 1 + X_target %*% w_true_target + rnorm(n_target, mean = 0, sd = sigma) ; Y_target = scale(Y_target, center = T, scale = F) ; Y_target_center = attr(Y_target, "scaled:center")
# 
# X_source = matrix(rnorm(n_source * p, mean = 0, sd = 1), n_source, p) ; X_source = scale(X_source, center = T, scale = F) ; X_source_center = attr(X_source, "scaled:center")
# Y_source = 1 + X_source %*% w_true_source + rnorm(n_source, sd = sigma); Y_source = scale(Y_source, center = T, scale = F) ; Y_source_center = attr(Y_source, "scaled:center")
# 
# 
# # cvxr
# cvxstl = cvx_stl(X0 = X_target, y0 = Y_target, X1 = X_source, y1 = Y_source,
#                   lambda0 = lambda0, lambda1 = lambda1)
# as.numeric(round(cvxstl$B0,4)) ;  as.numeric(round(cvxstl$B1,4))
# 
# # admm
# admmstl = admm_stl(X0 = X_target, y0 = Y_target, X1 = X_source, y1 = Y_source,
#                    lambda0 = lambda0, lambda1 = lambda1, tol_prim = 1e-8, tol_dual = 1e-8, verbose = T)
# as.numeric(round(admmstl$B0,4)) ;  as.numeric(round(admmstl$B1,4))
# 
# 
# # compare obj value
# cvxstl$result$value
# admmstl$objvalue[sum(!is.na(admmstl$objvalue))]

# Additional functions from functions_stl_cluster.R

admm_stl_sixsource <- function(X0, y0, X1, y1, X2, y2, X3, y3, X4, y4, X5, y5, X6, y6,
                               X0tX0 = NULL, X0ty0 = NULL,
                               X1tX1 = NULL, X1ty1 = NULL,
                               X2tX2 = NULL, X2ty2 = NULL,
                               X3tX3 = NULL, X3ty3 = NULL,
                               X4tX4 = NULL, X4ty4 = NULL,
                               X5tX5 = NULL, X5ty5 = NULL,
                               X6tX6 = NULL, X6ty6 = NULL,
                               lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6,
                               rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1,
                               rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1,
                               mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                               max_iter = 5000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3) + nrow(X4) + nrow(X5) + nrow(X6)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho02value = rep(NA, max_iter)
  rho03value = rep(NA, max_iter)
  rho04value = rep(NA, max_iter)
  rho05value = rep(NA, max_iter)
  rho06value = rep(NA, max_iter)
  
  rho1value = rep(NA, max_iter)
  rho2value = rep(NA, max_iter)
  rho3value = rep(NA, max_iter)
  rho4value = rep(NA, max_iter)
  rho5value = rep(NA, max_iter)
  rho6value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
  r_primal04value = rep(NA, max_iter) ; s_dual04value = rep(NA, max_iter)
  r_primal05value = rep(NA, max_iter) ; s_dual05value = rep(NA, max_iter)
  r_primal06value = rep(NA, max_iter) ; s_dual06value = rep(NA, max_iter)
  
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
  r_primal4value = rep(NA, max_iter) ; s_dual4value = rep(NA, max_iter)
  r_primal5value = rep(NA, max_iter) ; s_dual5value = rep(NA, max_iter)
  r_primal6value = rep(NA, max_iter) ; s_dual6value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  B2 <- matrix(0, p, 1)
  B3 <- matrix(0, p, 1)
  B4 <- matrix(0, p, 1)
  B5 <- matrix(0, p, 1)
  B6 <- matrix(0, p, 1)
  
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma02 <- matrix(0, p, 1)
  Gamma03 <- matrix(0, p, 1)
  Gamma04 <- matrix(0, p, 1)
  Gamma05 <- matrix(0, p, 1)
  Gamma06 <- matrix(0, p, 1)
  
  Gamma1 <- matrix(0, p, 1)
  Gamma2 <- matrix(0, p, 1)
  Gamma3 <- matrix(0, p, 1)
  Gamma4 <- matrix(0, p, 1)
  Gamma5 <- matrix(0, p, 1)
  Gamma6 <- matrix(0, p, 1)
  
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  U2 <- matrix(0, p, 1)
  U3 <- matrix(0, p, 1)
  U4 <- matrix(0, p, 1)
  U5 <- matrix(0, p, 1)
  U6 <- matrix(0, p, 1)
  
  V1 <- matrix(0, p, 1)
  V2 <- matrix(0, p, 1)
  V3 <- matrix(0, p, 1)
  V4 <- matrix(0, p, 1)
  V5 <- matrix(0, p, 1)
  V6 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0)
    X1tX1 = crossprod(X1,X1)
    X2tX2 = crossprod(X2,X2)
    X3tX3 = crossprod(X3,X3)
    X4tX4 = crossprod(X4,X4)
    X5tX5 = crossprod(X5,X5)
    X6tX6 = crossprod(X6,X6)
    X0ty0 = crossprod(X0,y0)
    X1ty1 = crossprod(X1,y1)
    X2ty2 = crossprod(X2,y2)
    X3ty3 = crossprod(X3,y3)
    X4ty4 = crossprod(X4,y4)
    X5ty5 = crossprod(X5,y5)
    X6ty6 = crossprod(X6,y6)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma03_old = Gamma03
    Gamma04_old = Gamma04
    Gamma05_old = Gamma05
    Gamma06_old = Gamma06
    
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    Gamma3_old = Gamma3
    Gamma4_old = Gamma4
    Gamma5_old = Gamma5
    Gamma6_old = Gamma6
    
    #### Step 1: B0 update ####
    lhs = 1/N*X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/N*X2tX2 + rho2*Ip
    rhs = 1/N*X2ty2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/N*X3tX3 + rho3*Ip
    rhs = 1/N*X3ty3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/N*X4tX4 + rho4*Ip
    rhs = 1/N*X4ty4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### B5 update ####
    lhs = 1/N*X5tX5 + rho5*Ip
    rhs = 1/N*X5ty5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    #### B6 update ####
    lhs = 1/N*X6tX6 + rho6*Ip
    rhs = 1/N*X6ty6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    
    #### Gamma00 update ####
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Gamma01 update ####
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    #### Gamma02 update ####
    Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
    #### Gamma03 update ####
    Gamma03 = Gamma3 + soft_threshold(B0+U3-Gamma3, lambda3/rho03)
    #### Gamma04 update ####
    Gamma04 = Gamma4 + soft_threshold(B0+U4-Gamma4, lambda4/rho04)
    #### Gamma05 update ####
    Gamma05 = Gamma5 + soft_threshold(B0+U5-Gamma5, lambda5/rho05)
    #### Gamma06 update ####
    Gamma06 = Gamma6 + soft_threshold(B0+U6-Gamma6, lambda6/rho06)
    
    #### Gamma1 update ####
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
    #### Gamma2 update ####
    Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
    #### Gamma3 update ####
    Gamma3 = Gamma03 + soft_threshold(B3+V3-Gamma03, lambda3/rho3)
    #### Gamma4 update ####
    Gamma4 = Gamma04 + soft_threshold(B4+V4-Gamma04, lambda4/rho4)
    #### Gamma5 update ####
    Gamma5 = Gamma05 + soft_threshold(B5+V5-Gamma05, lambda5/rho5)
    #### Gamma6 update ####
    Gamma6 = Gamma06 + soft_threshold(B6+V6-Gamma06, lambda6/rho6)
    
    #### Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    U3 <- U3 + B0 - Gamma03
    U4 <- U4 + B0 - Gamma04
    U5 <- U5 + B0 - Gamma05
    U6 <- U6 + B0 - Gamma06
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    V3 <- V3 + B3 - Gamma3
    V4 <- V4 + B4 - Gamma4
    V5 <- V5 + B5 - Gamma5
    V6 <- V6 + B6 - Gamma6
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
    r_primal03 = l2norm(B0 - Gamma03) ; s_dual03 = rho03 * l2norm(Gamma03 - Gamma03_old)
    r_primal04 = l2norm(B0 - Gamma04) ; s_dual04 = rho04 * l2norm(Gamma04 - Gamma04_old)
    r_primal05 = l2norm(B0 - Gamma05) ; s_dual05 = rho05 * l2norm(Gamma05 - Gamma05_old)
    r_primal06 = l2norm(B0 - Gamma06) ; s_dual06 = rho06 * l2norm(Gamma06 - Gamma06_old)
    
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
    r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
    r_primal3 = l2norm(B3 - Gamma3) ; s_dual3 = rho3 * l2norm(Gamma3 - Gamma3_old)
    r_primal4 = l2norm(B4 - Gamma4) ; s_dual4 = rho4 * l2norm(Gamma4 - Gamma4_old)
    r_primal5 = l2norm(B5 - Gamma5) ; s_dual5 = rho5 * l2norm(Gamma5 - Gamma5_old)
    r_primal6 = l2norm(B6 - Gamma6) ; s_dual6 = rho6 * l2norm(Gamma6 - Gamma6_old)
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
    r_primal03value[iter] = r_primal03 ; s_dual03value[iter] = s_dual03
    r_primal04value[iter] = r_primal04 ; s_dual04value[iter] = s_dual04
    r_primal05value[iter] = r_primal05 ; s_dual05value[iter] = s_dual05
    r_primal06value[iter] = r_primal06 ; s_dual06value[iter] = s_dual06
    
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
    r_primal3value[iter] = r_primal3 ; s_dual3value[iter] = s_dual3
    r_primal4value[iter] = r_primal4 ; s_dual4value[iter] = s_dual4
    r_primal5value[iter] = r_primal5 ; s_dual5value[iter] = s_dual5
    r_primal6value[iter] = r_primal6 ; s_dual6value[iter] = s_dual6
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal02 > mu * s_dual02){
      rho02 = min(rho02 * tau_incr, rho_max)
      U2 = U2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual02 > mu * r_primal02){
      rho02 = rho02 / tau_decr
      U2 = U2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal03 > mu * s_dual03){
      rho03 = min(rho03 * tau_incr, rho_max)
      U3 = U3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual03 > mu * r_primal03){
      rho03 = rho03 / tau_decr
      U3 = U3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal04 > mu * s_dual04){
      rho04 = min(rho04 * tau_incr, rho_max)
      U4 = U4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual04 > mu * r_primal04){
      rho04 = rho04 / tau_decr
      U4 = U4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal05 > mu * s_dual05){
      rho05 = min(rho05 * tau_incr, rho_max)
      U5 = U5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual05 > mu * r_primal05){
      rho05 = rho05 / tau_decr
      U5 = U5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal06 > mu * s_dual06){
      rho06 = min(rho06 * tau_incr, rho_max)
      U6 = U6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual06 > mu * r_primal06){
      rho06 = rho06 / tau_decr
      U6 = U6 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal2 > mu * s_dual2){
      rho2 = min(rho2 * tau_incr, rho_max)
      V2 = V2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual2 > mu * r_primal2){
      rho2 = rho2 / tau_decr
      V2 = V2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal3 > mu * s_dual3){
      rho3 = min(rho3 * tau_incr, rho_max)
      V3 = V3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual3 > mu * r_primal3){
      rho3 = rho3 / tau_decr
      V3 = V3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal4 > mu * s_dual4){
      rho4 = min(rho4 * tau_incr, rho_max)
      V4 = V4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual4 > mu * r_primal4){
      rho4 = rho4 / tau_decr
      V4 = V4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal5 > mu * s_dual5){
      rho5 = min(rho5 * tau_incr, rho_max)
      V5 = V5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual5 > mu * r_primal5){
      rho5 = rho5 / tau_decr
      V5 = V5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal6 > mu * s_dual6){
      rho6 = min(rho6 * tau_incr, rho_max)
      V6 = V6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual6 > mu * r_primal6){
      rho6 = rho6 / tau_decr
      V6 = V6 * tau_decr # change of variable required for the scaled variable
    }
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho02value[iter] = rho02
    rho03value[iter] = rho03
    rho04value[iter] = rho04
    rho05value[iter] = rho05
    rho06value[iter] = rho06
    
    rho1value[iter] = rho1
    rho2value[iter] = rho2
    rho3value[iter] = rho3
    rho4value[iter] = rho4
    rho5value[iter] = rho5
    rho6value[iter] = rho6
    
    if (max(r_primal00,
            r_primal01, r_primal02, r_primal03, r_primal04, r_primal05, r_primal06,
            r_primal1, r_primal2, r_primal3, r_primal4, r_primal5, r_primal6) < tol_prim
        &&
        max(s_dual00,
            s_dual01, s_dual02, s_dual03, s_dual04, s_dual05, s_dual06,
            s_dual1, s_dual2, s_dual3, s_dual4, s_dual5, s_dual6) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6,
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, rho04 = rho04value, rho05 = rho05value, rho06 = rho06value,
              rho1 = rho1value, rho2 = rho2value, rho3 = rho3value, rho4 = rho4value, rho5 = rho5value, rho6 = rho6value,
              
              r_primal00 = r_primal00value, s_dual00 = s_dual00value,
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              r_primal04 = r_primal04value, s_dual04 = s_dual04value,
              r_primal05 = r_primal05value, s_dual05 = s_dual05value,
              r_primal06 = r_primal06value, s_dual06 = s_dual06value,
              
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value,
              r_primal4 = r_primal4value, s_dual4 = s_dual4value,
              r_primal5 = r_primal5value, s_dual5 = s_dual5value,
              r_primal6 = r_primal6value, s_dual6 = s_dual6value))
}

admm_stl_elevensource <- function(X0, y0, X1, y1, X2, y2, X3, y3, X4, y4, X5, y5, X6, y6, X7, y7, X8, y8, X9, y9, X10, y10, X11, y11,
                                  X0tX0 = NULL, X0ty0 = NULL,
                                  X1tX1 = NULL, X1ty1 = NULL,
                                  X2tX2 = NULL, X2ty2 = NULL,
                                  X3tX3 = NULL, X3ty3 = NULL,
                                  X4tX4 = NULL, X4ty4 = NULL,
                                  X5tX5 = NULL, X5ty5 = NULL,
                                  X6tX6 = NULL, X6ty6 = NULL,
                                  X7tX7 = NULL, X7ty7 = NULL,
                                  X8tX8 = NULL, X8ty8 = NULL,
                                  X9tX9 = NULL, X9ty9 = NULL,
                                  X10tX10 = NULL, X10ty10 = NULL,
                                  X11tX11 = NULL, X11ty11 = NULL,
                                  lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8, lambda9, lambda10, lambda11,
                                  rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1, rho07 = 1, rho08 = 1, rho09 = 1, rho010 = 1, rho011 = 1,
                                  rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1, rho7 = 1, rho8 = 1, rho9 = 1, rho10 = 1, rho11 = 1,
                                  mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                  max_iter = 5000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3) + nrow(X4) + nrow(X5) + nrow(X6) + nrow(X7) + nrow(X8) + nrow(X9) + nrow(X10) + nrow(X11)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho02value = rep(NA, max_iter)
  rho03value = rep(NA, max_iter)
  rho04value = rep(NA, max_iter)
  rho05value = rep(NA, max_iter)
  rho06value = rep(NA, max_iter)
  rho07value = rep(NA, max_iter)
  rho08value = rep(NA, max_iter)
  rho09value = rep(NA, max_iter)
  rho010value = rep(NA, max_iter)
  rho011value = rep(NA, max_iter)
  
  rho1value = rep(NA, max_iter)
  rho2value = rep(NA, max_iter)
  rho3value = rep(NA, max_iter)
  rho4value = rep(NA, max_iter)
  rho5value = rep(NA, max_iter)
  rho6value = rep(NA, max_iter)
  rho7value = rep(NA, max_iter)
  rho8value = rep(NA, max_iter)
  rho9value = rep(NA, max_iter)
  rho10value = rep(NA, max_iter)
  rho11value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
  r_primal04value = rep(NA, max_iter) ; s_dual04value = rep(NA, max_iter)
  r_primal05value = rep(NA, max_iter) ; s_dual05value = rep(NA, max_iter)
  r_primal06value = rep(NA, max_iter) ; s_dual06value = rep(NA, max_iter)
  r_primal07value = rep(NA, max_iter) ; s_dual07value = rep(NA, max_iter)
  r_primal08value = rep(NA, max_iter) ; s_dual08value = rep(NA, max_iter)
  r_primal09value = rep(NA, max_iter) ; s_dual09value = rep(NA, max_iter)
  r_primal010value = rep(NA, max_iter) ; s_dual010value = rep(NA, max_iter)
  r_primal011value = rep(NA, max_iter) ; s_dual011value = rep(NA, max_iter)
  
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
  r_primal4value = rep(NA, max_iter) ; s_dual4value = rep(NA, max_iter)
  r_primal5value = rep(NA, max_iter) ; s_dual5value = rep(NA, max_iter)
  r_primal6value = rep(NA, max_iter) ; s_dual6value = rep(NA, max_iter)
  r_primal7value = rep(NA, max_iter) ; s_dual7value = rep(NA, max_iter)
  r_primal8value = rep(NA, max_iter) ; s_dual8value = rep(NA, max_iter)
  r_primal9value = rep(NA, max_iter) ; s_dual9value = rep(NA, max_iter)
  r_primal10value = rep(NA, max_iter) ; s_dual10value = rep(NA, max_iter)
  r_primal11value = rep(NA, max_iter) ; s_dual11value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  B2 <- matrix(0, p, 1)
  B3 <- matrix(0, p, 1)
  B4 <- matrix(0, p, 1)
  B5 <- matrix(0, p, 1)
  B6 <- matrix(0, p, 1)
  B7 <- matrix(0, p, 1)
  B8 <- matrix(0, p, 1)
  B9 <- matrix(0, p, 1)
  B10 <- matrix(0, p, 1)
  B11 <- matrix(0, p, 1)
  
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma02 <- matrix(0, p, 1)
  Gamma03 <- matrix(0, p, 1)
  Gamma04 <- matrix(0, p, 1)
  Gamma05 <- matrix(0, p, 1)
  Gamma06 <- matrix(0, p, 1)
  Gamma07 <- matrix(0, p, 1)
  Gamma08 <- matrix(0, p, 1)
  Gamma09 <- matrix(0, p, 1)
  Gamma010 <- matrix(0, p, 1)
  Gamma011 <- matrix(0, p, 1)
  
  Gamma1 <- matrix(0, p, 1)
  Gamma2 <- matrix(0, p, 1)
  Gamma3 <- matrix(0, p, 1)
  Gamma4 <- matrix(0, p, 1)
  Gamma5 <- matrix(0, p, 1)
  Gamma6 <- matrix(0, p, 1)
  Gamma7 <- matrix(0, p, 1)
  Gamma8 <- matrix(0, p, 1)
  Gamma9 <- matrix(0, p, 1)
  Gamma10 <- matrix(0, p, 1)
  Gamma11 <- matrix(0, p, 1)
  
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  U2 <- matrix(0, p, 1)
  U3 <- matrix(0, p, 1)
  U4 <- matrix(0, p, 1)
  U5 <- matrix(0, p, 1)
  U6 <- matrix(0, p, 1)
  U7 <- matrix(0, p, 1)
  U8 <- matrix(0, p, 1)
  U9 <- matrix(0, p, 1)
  U10 <- matrix(0, p, 1)
  U11 <- matrix(0, p, 1)
  
  V1 <- matrix(0, p, 1)
  V2 <- matrix(0, p, 1)
  V3 <- matrix(0, p, 1)
  V4 <- matrix(0, p, 1)
  V5 <- matrix(0, p, 1)
  V6 <- matrix(0, p, 1)
  V7 <- matrix(0, p, 1)
  V8 <- matrix(0, p, 1)
  V9 <- matrix(0, p, 1)
  V10 <- matrix(0, p, 1)
  V11 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0)
    X1tX1 = crossprod(X1,X1)
    X2tX2 = crossprod(X2,X2)
    X3tX3 = crossprod(X3,X3)
    X4tX4 = crossprod(X4,X4)
    X5tX5 = crossprod(X5,X5)
    X6tX6 = crossprod(X6,X6)
    X7tX7 = crossprod(X7,X7)
    X8tX8 = crossprod(X8,X8)
    X9tX9 = crossprod(X9,X9)
    X10tX10 = crossprod(X10,X10)
    X11tX11 = crossprod(X11,X11)
    X0ty0 = crossprod(X0,y0)
    X1ty1 = crossprod(X1,y1)
    X2ty2 = crossprod(X2,y2)
    X3ty3 = crossprod(X3,y3)
    X4ty4 = crossprod(X4,y4)
    X5ty5 = crossprod(X5,y5)
    X6ty6 = crossprod(X6,y6)
    X7ty7 = crossprod(X7,y7)
    X8ty8 = crossprod(X8,y8)
    X9ty9 = crossprod(X9,y9)
    X10ty10 = crossprod(X10,y10)
    X11ty11 = crossprod(X11,y11)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma03_old = Gamma03
    Gamma04_old = Gamma04
    Gamma05_old = Gamma05
    Gamma06_old = Gamma06
    Gamma07_old = Gamma07
    Gamma08_old = Gamma08
    Gamma09_old = Gamma09
    Gamma010_old = Gamma010
    Gamma011_old = Gamma011
    
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    Gamma3_old = Gamma3
    Gamma4_old = Gamma4
    Gamma5_old = Gamma5
    Gamma6_old = Gamma6
    Gamma7_old = Gamma7
    Gamma8_old = Gamma8
    Gamma9_old = Gamma9
    Gamma10_old = Gamma10
    Gamma11_old = Gamma11
    
    #### Step 1: B0 update ####
    lhs = 1/N*X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06 + rho07 + rho08 + rho09 + rho010 + rho011)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6) + rho07 * (Gamma07 - U7) + rho08 * (Gamma08 - U8) + rho09 * (Gamma09 - U9) + rho010 * (Gamma010 - U10) + rho011 * (Gamma011 - U11)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/N*X2tX2 + rho2*Ip
    rhs = 1/N*X2ty2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/N*X3tX3 + rho3*Ip
    rhs = 1/N*X3ty3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/N*X4tX4 + rho4*Ip
    rhs = 1/N*X4ty4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### B5 update ####
    lhs = 1/N*X5tX5 + rho5*Ip
    rhs = 1/N*X5ty5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    #### B6 update ####
    lhs = 1/N*X6tX6 + rho6*Ip
    rhs = 1/N*X6ty6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    
    #### B7 update ####
    lhs = 1/N*X7tX7 + rho7*Ip
    rhs = 1/N*X7ty7 + rho7*(Gamma7 - V7)
    B7 <- solve(lhs, rhs)
    
    #### B8 update ####
    lhs = 1/N*X8tX8 + rho8*Ip
    rhs = 1/N*X8ty8 + rho8*(Gamma8 - V8)
    B8 <- solve(lhs, rhs)
    
    #### B9 update ####
    lhs = 1/N*X9tX9 + rho9*Ip
    rhs = 1/N*X9ty9 + rho9*(Gamma9 - V9)
    B9 <- solve(lhs, rhs)
    
    #### B10 update ####
    lhs = 1/N*X10tX10 + rho10*Ip
    rhs = 1/N*X10ty10 + rho10*(Gamma10 - V10)
    B10 <- solve(lhs, rhs)
    
    #### B11 update ####
    lhs = 1/N*X11tX11 + rho11*Ip
    rhs = 1/N*X11ty11 + rho11*(Gamma11 - V11)
    B11 <- solve(lhs, rhs)
    
    #### Gamma00 update ####
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Gamma01 update ####
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    #### Gamma02 update ####
    Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
    #### Gamma03 update ####
    Gamma03 = Gamma3 + soft_threshold(B0+U3-Gamma3, lambda3/rho03)
    #### Gamma04 update ####
    Gamma04 = Gamma4 + soft_threshold(B0+U4-Gamma4, lambda4/rho04)
    #### Gamma05 update ####
    Gamma05 = Gamma5 + soft_threshold(B0+U5-Gamma5, lambda5/rho05)
    #### Gamma06 update ####
    Gamma06 = Gamma6 + soft_threshold(B0+U6-Gamma6, lambda6/rho06)
    #### Gamma07 update ####
    Gamma07 = Gamma7 + soft_threshold(B0+U7-Gamma7, lambda7/rho07)
    #### Gamma08 update ####
    Gamma08 = Gamma8 + soft_threshold(B0+U8-Gamma8, lambda8/rho08)
    #### Gamma09 update ####
    Gamma09 = Gamma9 + soft_threshold(B0+U9-Gamma9, lambda9/rho09)
    #### Gamma010 update ####
    Gamma010 = Gamma10 + soft_threshold(B0+U10-Gamma10, lambda10/rho010)
    #### Gamma011 update ####
    Gamma011 = Gamma11 + soft_threshold(B0+U11-Gamma11, lambda11/rho011)
    
    #### Gamma1 update ####
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
    #### Gamma2 update ####
    Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
    #### Gamma3 update ####
    Gamma3 = Gamma03 + soft_threshold(B3+V3-Gamma03, lambda3/rho3)
    #### Gamma4 update ####
    Gamma4 = Gamma04 + soft_threshold(B4+V4-Gamma04, lambda4/rho4)
    #### Gamma5 update ####
    Gamma5 = Gamma05 + soft_threshold(B5+V5-Gamma05, lambda5/rho5)
    #### Gamma6 update ####
    Gamma6 = Gamma06 + soft_threshold(B6+V6-Gamma06, lambda6/rho6)
    #### Gamma7 update ####
    Gamma7 = Gamma07 + soft_threshold(B7+V7-Gamma07, lambda7/rho7)
    #### Gamma8 update ####
    Gamma8 = Gamma08 + soft_threshold(B8+V8-Gamma08, lambda8/rho8)
    #### Gamma9 update ####
    Gamma9 = Gamma09 + soft_threshold(B9+V9-Gamma09, lambda9/rho9)
    #### Gamma10 update ####
    Gamma10 = Gamma010 + soft_threshold(B10+V10-Gamma010, lambda10/rho10)
    #### Gamma11 update ####
    Gamma11 = Gamma011 + soft_threshold(B11+V11-Gamma011, lambda11/rho11)
    
    #### Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    U3 <- U3 + B0 - Gamma03
    U4 <- U4 + B0 - Gamma04
    U5 <- U5 + B0 - Gamma05
    U6 <- U6 + B0 - Gamma06
    U7 <- U7 + B0 - Gamma07
    U8 <- U8 + B0 - Gamma08
    U9 <- U9 + B0 - Gamma09
    U10 <- U10 + B0 - Gamma010
    U11 <- U11 + B0 - Gamma011
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    V3 <- V3 + B3 - Gamma3
    V4 <- V4 + B4 - Gamma4
    V5 <- V5 + B5 - Gamma5
    V6 <- V6 + B6 - Gamma6
    V7 <- V7 + B7 - Gamma7
    V8 <- V8 + B8 - Gamma8
    V9 <- V9 + B9 - Gamma9
    V10 <- V10 + B10 - Gamma10
    V11 <- V11 + B11 - Gamma11
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
    r_primal03 = l2norm(B0 - Gamma03) ; s_dual03 = rho03 * l2norm(Gamma03 - Gamma03_old)
    r_primal04 = l2norm(B0 - Gamma04) ; s_dual04 = rho04 * l2norm(Gamma04 - Gamma04_old)
    r_primal05 = l2norm(B0 - Gamma05) ; s_dual05 = rho05 * l2norm(Gamma05 - Gamma05_old)
    r_primal06 = l2norm(B0 - Gamma06) ; s_dual06 = rho06 * l2norm(Gamma06 - Gamma06_old)
    r_primal07 = l2norm(B0 - Gamma07) ; s_dual07 = rho07 * l2norm(Gamma07 - Gamma07_old)
    r_primal08 = l2norm(B0 - Gamma08) ; s_dual08 = rho08 * l2norm(Gamma08 - Gamma08_old)
    r_primal09 = l2norm(B0 - Gamma09) ; s_dual09 = rho09 * l2norm(Gamma09 - Gamma09_old)
    r_primal010 = l2norm(B0 - Gamma010) ; s_dual010 = rho010 * l2norm(Gamma010 - Gamma010_old)
    r_primal011 = l2norm(B0 - Gamma011) ; s_dual011 = rho011 * l2norm(Gamma011 - Gamma011_old)
    
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
    r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
    r_primal3 = l2norm(B3 - Gamma3) ; s_dual3 = rho3 * l2norm(Gamma3 - Gamma3_old)
    r_primal4 = l2norm(B4 - Gamma4) ; s_dual4 = rho4 * l2norm(Gamma4 - Gamma4_old)
    r_primal5 = l2norm(B5 - Gamma5) ; s_dual5 = rho5 * l2norm(Gamma5 - Gamma5_old)
    r_primal6 = l2norm(B6 - Gamma6) ; s_dual6 = rho6 * l2norm(Gamma6 - Gamma6_old)
    r_primal7 = l2norm(B7 - Gamma7) ; s_dual7 = rho7 * l2norm(Gamma7 - Gamma7_old)
    r_primal8 = l2norm(B8 - Gamma8) ; s_dual8 = rho8 * l2norm(Gamma8 - Gamma8_old)
    r_primal9 = l2norm(B9 - Gamma9) ; s_dual9 = rho9 * l2norm(Gamma9 - Gamma9_old)
    r_primal10 = l2norm(B10 - Gamma10) ; s_dual10 = rho10 * l2norm(Gamma10 - Gamma10_old)
    r_primal11 = l2norm(B11 - Gamma11) ; s_dual11 = rho11 * l2norm(Gamma11 - Gamma11_old)
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
    r_primal03value[iter] = r_primal03 ; s_dual03value[iter] = s_dual03
    r_primal04value[iter] = r_primal04 ; s_dual04value[iter] = s_dual04
    r_primal05value[iter] = r_primal05 ; s_dual05value[iter] = s_dual05
    r_primal06value[iter] = r_primal06 ; s_dual06value[iter] = s_dual06
    r_primal07value[iter] = r_primal07 ; s_dual07value[iter] = s_dual07
    r_primal08value[iter] = r_primal08 ; s_dual08value[iter] = s_dual08
    r_primal09value[iter] = r_primal09 ; s_dual09value[iter] = s_dual09
    r_primal010value[iter] = r_primal010 ; s_dual010value[iter] = s_dual010
    r_primal011value[iter] = r_primal011 ; s_dual011value[iter] = s_dual011
    
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
    r_primal3value[iter] = r_primal3 ; s_dual3value[iter] = s_dual3
    r_primal4value[iter] = r_primal4 ; s_dual4value[iter] = s_dual4
    r_primal5value[iter] = r_primal5 ; s_dual5value[iter] = s_dual5
    r_primal6value[iter] = r_primal6 ; s_dual6value[iter] = s_dual6
    r_primal7value[iter] = r_primal7 ; s_dual7value[iter] = s_dual7
    r_primal8value[iter] = r_primal8 ; s_dual8value[iter] = s_dual8
    r_primal9value[iter] = r_primal9 ; s_dual9value[iter] = s_dual9
    r_primal10value[iter] = r_primal10 ; s_dual10value[iter] = s_dual10
    r_primal11value[iter] = r_primal11 ; s_dual11value[iter] = s_dual11
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal02 > mu * s_dual02){
      rho02 = min(rho02 * tau_incr, rho_max)
      U2 = U2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual02 > mu * r_primal02){
      rho02 = rho02 / tau_decr
      U2 = U2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal03 > mu * s_dual03){
      rho03 = min(rho03 * tau_incr, rho_max)
      U3 = U3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual03 > mu * r_primal03){
      rho03 = rho03 / tau_decr
      U3 = U3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal04 > mu * s_dual04){
      rho04 = min(rho04 * tau_incr, rho_max)
      U4 = U4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual04 > mu * r_primal04){
      rho04 = rho04 / tau_decr
      U4 = U4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal05 > mu * s_dual05){
      rho05 = min(rho05 * tau_incr, rho_max)
      U5 = U5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual05 > mu * r_primal05){
      rho05 = rho05 / tau_decr
      U5 = U5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal06 > mu * s_dual06){
      rho06 = min(rho06 * tau_incr, rho_max)
      U6 = U6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual06 > mu * r_primal06){
      rho06 = rho06 / tau_decr
      U6 = U6 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal07 > mu * s_dual07){
      rho07 = min(rho07 * tau_incr, rho_max)
      U7 = U7 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual07 > mu * r_primal07){
      rho07 = rho07 / tau_decr
      U7 = U7 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal08 > mu * s_dual08){
      rho08 = min(rho08 * tau_incr, rho_max)
      U8 = U8 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual08 > mu * r_primal08){
      rho08 = rho08 / tau_decr
      U8 = U8 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal09 > mu * s_dual09){
      rho09 = min(rho09 * tau_incr, rho_max)
      U9 = U9 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual09 > mu * r_primal09){
      rho09 = rho09 / tau_decr
      U9 = U9 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal010 > mu * s_dual010){
      rho010 = min(rho010 * tau_incr, rho_max)
      U10 = U10 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual010 > mu * r_primal010){
      rho010 = rho010 / tau_decr
      U10 = U10 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal011 > mu * s_dual011){
      rho011 = min(rho011 * tau_incr, rho_max)
      U11 = U11 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual011 > mu * r_primal011){
      rho011 = rho011 / tau_decr
      U11 = U11 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal2 > mu * s_dual2){
      rho2 = min(rho2 * tau_incr, rho_max)
      V2 = V2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual2 > mu * r_primal2){
      rho2 = rho2 / tau_decr
      V2 = V2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal3 > mu * s_dual3){
      rho3 = min(rho3 * tau_incr, rho_max)
      V3 = V3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual3 > mu * r_primal3){
      rho3 = rho3 / tau_decr
      V3 = V3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal4 > mu * s_dual4){
      rho4 = min(rho4 * tau_incr, rho_max)
      V4 = V4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual4 > mu * r_primal4){
      rho4 = rho4 / tau_decr
      V4 = V4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal5 > mu * s_dual5){
      rho5 = min(rho5 * tau_incr, rho_max)
      V5 = V5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual5 > mu * r_primal5){
      rho5 = rho5 / tau_decr
      V5 = V5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal6 > mu * s_dual6){
      rho6 = min(rho6 * tau_incr, rho_max)
      V6 = V6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual6 > mu * r_primal6){
      rho6 = rho6 / tau_decr
      V6 = V6 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal7 > mu * s_dual7){
      rho7 = min(rho7 * tau_incr, rho_max)
      V7 = V7 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual7 > mu * r_primal7){
      rho7 = rho7 / tau_decr
      V7 = V7 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal8 > mu * s_dual8){
      rho8 = min(rho8 * tau_incr, rho_max)
      V8 = V8 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual8 > mu * r_primal8){
      rho8 = rho8 / tau_decr
      V8 = V8 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal9 > mu * s_dual9){
      rho9 = min(rho9 * tau_incr, rho_max)
      V9 = V9 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual9 > mu * r_primal9){
      rho9 = rho9 / tau_decr
      V9 = V9 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal10 > mu * s_dual10){
      rho10 = min(rho10 * tau_incr, rho_max)
      V10 = V10 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual10 > mu * r_primal10){
      rho10 = rho10 / tau_decr
      V10 = V10 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal11 > mu * s_dual11){
      rho11 = min(rho11 * tau_incr, rho_max)
      V11 = V11 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual11 > mu * r_primal11){
      rho11 = rho11 / tau_decr
      V11 = V11 * tau_decr # change of variable required for the scaled variable
    }
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho02value[iter] = rho02
    rho03value[iter] = rho03
    rho04value[iter] = rho04
    rho05value[iter] = rho05
    rho06value[iter] = rho06
    rho07value[iter] = rho07
    rho08value[iter] = rho08
    rho09value[iter] = rho09
    rho010value[iter] = rho010
    rho011value[iter] = rho011
    
    rho1value[iter] = rho1
    rho2value[iter] = rho2
    rho3value[iter] = rho3
    rho4value[iter] = rho4
    rho5value[iter] = rho5
    rho6value[iter] = rho6
    rho7value[iter] = rho7
    rho8value[iter] = rho8
    rho9value[iter] = rho9
    rho10value[iter] = rho10
    rho11value[iter] = rho11
    
    if (max(r_primal00,
            r_primal01, r_primal02, r_primal03, r_primal04, r_primal05, r_primal06, r_primal07, r_primal08, r_primal09, r_primal010, r_primal011,
            r_primal1, r_primal2, r_primal3, r_primal4, r_primal5, r_primal6, r_primal7, r_primal8, r_primal9, r_primal10, r_primal11) < tol_prim
        &&
        max(s_dual00,
            s_dual01, s_dual02, s_dual03, s_dual04, s_dual05, s_dual06, s_dual07, s_dual08, s_dual09, s_dual010, s_dual011,
            s_dual1, s_dual2, s_dual3, s_dual4, s_dual5, s_dual6, s_dual7, s_dual8, s_dual9, s_dual10, s_dual11) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6, B7 = B7, B8 = B8, B9 = B9, B10 = B10, B11 = B11,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06, Gamma07 = Gamma07, Gamma08 = Gamma08, Gamma09 = Gamma09, Gamma010 = Gamma010, Gamma011 = Gamma011,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6, Gamma7 = Gamma7, Gamma8 = Gamma8, Gamma9 = Gamma9, Gamma10 = Gamma10, Gamma11 = Gamma11,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6, U7 = U7, U8 = U8, U9 = U9, U10 = U10, U11 = U11,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6, V7 = V7, V8 = V8, V9 = V9, V10 = V10, V11 = V11,
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, rho04 = rho04value, rho05 = rho05value, rho06 = rho06value, rho07 = rho07value, rho08 = rho08value, rho09 = rho09value, rho010 = rho010value, rho011 = rho011value,
              rho1 = rho1value, rho2 = rho2value, rho3 = rho3value, rho4 = rho4value, rho5 = rho5value, rho6 = rho6value, rho7 = rho7value, rho8 = rho8value, rho9 = rho9value, rho10 = rho10value, rho11 = rho11value,
              
              r_primal00 = r_primal00value, s_dual00 = s_dual00value,
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              r_primal04 = r_primal04value, s_dual04 = s_dual04value,
              r_primal05 = r_primal05value, s_dual05 = s_dual05value,
              r_primal06 = r_primal06value, s_dual06 = s_dual06value,
              r_primal07 = r_primal07value, s_dual07 = s_dual07value,
              r_primal08 = r_primal08value, s_dual08 = s_dual08value,
              r_primal09 = r_primal09value, s_dual09 = s_dual09value,
              r_primal010 = r_primal010value, s_dual010 = s_dual010value,
              r_primal011 = r_primal011value, s_dual011 = s_dual011value,
              
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value,
              r_primal4 = r_primal4value, s_dual4 = s_dual4value,
              r_primal5 = r_primal5value, s_dual5 = s_dual5value,
              r_primal6 = r_primal6value, s_dual6 = s_dual6value,
              r_primal7 = r_primal7value, s_dual7 = s_dual7value,
              r_primal8 = r_primal8value, s_dual8 = s_dual8value,
              r_primal9 = r_primal9value, s_dual9 = s_dual9value,
              r_primal10 = r_primal10value, s_dual10 = s_dual10value,
              r_primal11 = r_primal11value, s_dual11 = s_dual11value))
}

admm_stl_foursource <- function(X0, y0, X1, y1, X2, y2, X3, y3, X4, y4,
                                X0tX0 = NULL, X0ty0 = NULL,
                                X1tX1 = NULL, X1ty1 = NULL,
                                X2tX2 = NULL, X2ty2 = NULL,
                                X3tX3 = NULL, X3ty3 = NULL,
                                X4tX4 = NULL, X4ty4 = NULL,
                                lambda0, lambda1, lambda2, lambda3, lambda4,
                                rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1,
                                rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1,
                                mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                max_iter = 5000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3) + nrow(X4)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho02value = rep(NA, max_iter)
  rho03value = rep(NA, max_iter)
  rho04value = rep(NA, max_iter)
  
  rho1value = rep(NA, max_iter)
  rho2value = rep(NA, max_iter)
  rho3value = rep(NA, max_iter)
  rho4value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
  r_primal04value = rep(NA, max_iter) ; s_dual04value = rep(NA, max_iter)
  
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
  r_primal4value = rep(NA, max_iter) ; s_dual4value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  B2 <- matrix(0, p, 1)
  B3 <- matrix(0, p, 1)
  B4 <- matrix(0, p, 1)
  
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma02 <- matrix(0, p, 1)
  Gamma03 <- matrix(0, p, 1)
  Gamma04 <- matrix(0, p, 1)
  
  Gamma1 <- matrix(0, p, 1)
  Gamma2 <- matrix(0, p, 1)
  Gamma3 <- matrix(0, p, 1)
  Gamma4 <- matrix(0, p, 1)
  
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  U2 <- matrix(0, p, 1)
  U3 <- matrix(0, p, 1)
  U4 <- matrix(0, p, 1)
  
  V1 <- matrix(0, p, 1)
  V2 <- matrix(0, p, 1)
  V3 <- matrix(0, p, 1)
  V4 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0)
    X1tX1 = crossprod(X1,X1)
    X2tX2 = crossprod(X2,X2)
    X3tX3 = crossprod(X3,X3)
    X4tX4 = crossprod(X4,X4)
    X0ty0 = crossprod(X0,y0)
    X1ty1 = crossprod(X1,y1)
    X2ty2 = crossprod(X2,y2)
    X3ty3 = crossprod(X3,y3)
    X4ty4 = crossprod(X4,y4)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma03_old = Gamma03
    Gamma04_old = Gamma04
    
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    Gamma3_old = Gamma3
    Gamma4_old = Gamma4
    
    #### Step 1: B0 update ####
    lhs = 1/N*X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/N*X2tX2 + rho2*Ip
    rhs = 1/N*X2ty2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/N*X3tX3 + rho3*Ip
    rhs = 1/N*X3ty3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/N*X4tX4 + rho4*Ip
    rhs = 1/N*X4ty4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### Gamma00 update ####
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Gamma01 update ####
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    #### Gamma02 update ####
    Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
    #### Gamma03 update ####
    Gamma03 = Gamma3 + soft_threshold(B0+U3-Gamma3, lambda3/rho03)
    #### Gamma04 update ####
    Gamma04 = Gamma4 + soft_threshold(B0+U4-Gamma4, lambda4/rho04)
    
    #### Gamma1 update ####
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
    #### Gamma2 update ####
    Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
    #### Gamma3 update ####
    Gamma3 = Gamma03 + soft_threshold(B3+V3-Gamma03, lambda3/rho3)
    #### Gamma4 update ####
    Gamma4 = Gamma04 + soft_threshold(B4+V4-Gamma04, lambda4/rho4)
    
    #### Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    U3 <- U3 + B0 - Gamma03
    U4 <- U4 + B0 - Gamma04
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    V3 <- V3 + B3 - Gamma3
    V4 <- V4 + B4 - Gamma4
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
    r_primal03 = l2norm(B0 - Gamma03) ; s_dual03 = rho03 * l2norm(Gamma03 - Gamma03_old)
    r_primal04 = l2norm(B0 - Gamma04) ; s_dual04 = rho04 * l2norm(Gamma04 - Gamma04_old)
    
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
    r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
    r_primal3 = l2norm(B3 - Gamma3) ; s_dual3 = rho3 * l2norm(Gamma3 - Gamma3_old)
    r_primal4 = l2norm(B4 - Gamma4) ; s_dual4 = rho4 * l2norm(Gamma4 - Gamma4_old)
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
    r_primal03value[iter] = r_primal03 ; s_dual03value[iter] = s_dual03
    r_primal04value[iter] = r_primal04 ; s_dual04value[iter] = s_dual04
    
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
    r_primal3value[iter] = r_primal3 ; s_dual3value[iter] = s_dual3
    r_primal4value[iter] = r_primal4 ; s_dual4value[iter] = s_dual4
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal02 > mu * s_dual02){
      rho02 = min(rho02 * tau_incr, rho_max)
      U2 = U2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual02 > mu * r_primal02){
      rho02 = rho02 / tau_decr
      U2 = U2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal03 > mu * s_dual03){
      rho03 = min(rho03 * tau_incr, rho_max)
      U3 = U3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual03 > mu * r_primal03){
      rho03 = rho03 / tau_decr
      U3 = U3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal04 > mu * s_dual04){
      rho04 = min(rho04 * tau_incr, rho_max)
      U4 = U4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual04 > mu * r_primal04){
      rho04 = rho04 / tau_decr
      U4 = U4 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal2 > mu * s_dual2){
      rho2 = min(rho2 * tau_incr, rho_max)
      V2 = V2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual2 > mu * r_primal2){
      rho2 = rho2 / tau_decr
      V2 = V2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal3 > mu * s_dual3){
      rho3 = min(rho3 * tau_incr, rho_max)
      V3 = V3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual3 > mu * r_primal3){
      rho3 = rho3 / tau_decr
      V3 = V3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal4 > mu * s_dual4){
      rho4 = min(rho4 * tau_incr, rho_max)
      V4 = V4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual4 > mu * r_primal4){
      rho4 = rho4 / tau_decr
      V4 = V4 * tau_decr # change of variable required for the scaled variable
    }
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho02value[iter] = rho02
    rho03value[iter] = rho03
    rho04value[iter] = rho04
    
    rho1value[iter] = rho1
    rho2value[iter] = rho2
    rho3value[iter] = rho3
    rho4value[iter] = rho4
    
    if (max(r_primal00,
            r_primal01, r_primal02, r_primal03, r_primal04,
            r_primal1, r_primal2, r_primal3, r_primal4) < tol_prim
        &&
        max(s_dual00,
            s_dual01, s_dual02, s_dual03, s_dual04,
            s_dual1, s_dual2, s_dual3, s_dual4) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4,
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, rho04 = rho04value,
              rho1 = rho1value, rho2 = rho2value, rho3 = rho3value, rho4 = rho4value,
              
              r_primal00 = r_primal00value, s_dual00 = s_dual00value,
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              r_primal04 = r_primal04value, s_dual04 = s_dual04value,
              
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value,
              r_primal4 = r_primal4value, s_dual4 = s_dual4value))
}

admm_stl_thirteensource <- function(X0, y0, X1, y1, X2, y2, X3, y3, X4, y4, X5, y5, X6, y6, X7, y7, X8, y8, X9, y9, X10, y10, X11, y11, X12, y12, X13, y13,
                                    X0tX0 = NULL, X0ty0 = NULL,
                                    X1tX1 = NULL, X1ty1 = NULL,
                                    X2tX2 = NULL, X2ty2 = NULL,
                                    X3tX3 = NULL, X3ty3 = NULL,
                                    X4tX4 = NULL, X4ty4 = NULL,
                                    X5tX5 = NULL, X5ty5 = NULL,
                                    X6tX6 = NULL, X6ty6 = NULL,
                                    X7tX7 = NULL, X7ty7 = NULL,
                                    X8tX8 = NULL, X8ty8 = NULL,
                                    X9tX9 = NULL, X9ty9 = NULL,
                                    X10tX10 = NULL, X10ty10 = NULL,
                                    X11tX11 = NULL, X11ty11 = NULL,
                                    X12tX12 = NULL, X12ty12 = NULL,
                                    X13tX13 = NULL, X13ty13 = NULL,
                                    lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8, lambda9, lambda10, lambda11, lambda12, lambda13,
                                    rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1, rho07 = 1, rho08 = 1, rho09 = 1, rho010 = 1, rho011 = 1, rho012 = 1, rho013 = 1,
                                    rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1, rho7 = 1, rho8 = 1, rho9 = 1, rho10 = 1, rho11 = 1, rho12 = 1, rho13 = 1,
                                    mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                    max_iter = 5000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3) + nrow(X4) + nrow(X5) + nrow(X6) + nrow(X7) + nrow(X8) + nrow(X9) + nrow(X10) + nrow(X11) + nrow(X12) + nrow(X13)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho02value = rep(NA, max_iter)
  rho03value = rep(NA, max_iter)
  rho04value = rep(NA, max_iter)
  rho05value = rep(NA, max_iter)
  rho06value = rep(NA, max_iter)
  rho07value = rep(NA, max_iter)
  rho08value = rep(NA, max_iter)
  rho09value = rep(NA, max_iter)
  rho010value = rep(NA, max_iter)
  rho011value = rep(NA, max_iter)
  rho012value = rep(NA, max_iter)
  rho013value = rep(NA, max_iter)
  
  rho1value = rep(NA, max_iter)
  rho2value = rep(NA, max_iter)
  rho3value = rep(NA, max_iter)
  rho4value = rep(NA, max_iter)
  rho5value = rep(NA, max_iter)
  rho6value = rep(NA, max_iter)
  rho7value = rep(NA, max_iter)
  rho8value = rep(NA, max_iter)
  rho9value = rep(NA, max_iter)
  rho10value = rep(NA, max_iter)
  rho11value = rep(NA, max_iter)
  rho12value = rep(NA, max_iter)
  rho13value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
  r_primal04value = rep(NA, max_iter) ; s_dual04value = rep(NA, max_iter)
  r_primal05value = rep(NA, max_iter) ; s_dual05value = rep(NA, max_iter)
  r_primal06value = rep(NA, max_iter) ; s_dual06value = rep(NA, max_iter)
  r_primal07value = rep(NA, max_iter) ; s_dual07value = rep(NA, max_iter)
  r_primal08value = rep(NA, max_iter) ; s_dual08value = rep(NA, max_iter)
  r_primal09value = rep(NA, max_iter) ; s_dual09value = rep(NA, max_iter)
  r_primal010value = rep(NA, max_iter) ; s_dual010value = rep(NA, max_iter)
  r_primal011value = rep(NA, max_iter) ; s_dual011value = rep(NA, max_iter)
  r_primal012value = rep(NA, max_iter) ; s_dual012value = rep(NA, max_iter)
  r_primal013value = rep(NA, max_iter) ; s_dual013value = rep(NA, max_iter)
  
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
  r_primal4value = rep(NA, max_iter) ; s_dual4value = rep(NA, max_iter)
  r_primal5value = rep(NA, max_iter) ; s_dual5value = rep(NA, max_iter)
  r_primal6value = rep(NA, max_iter) ; s_dual6value = rep(NA, max_iter)
  r_primal7value = rep(NA, max_iter) ; s_dual7value = rep(NA, max_iter)
  r_primal8value = rep(NA, max_iter) ; s_dual8value = rep(NA, max_iter)
  r_primal9value = rep(NA, max_iter) ; s_dual9value = rep(NA, max_iter)
  r_primal10value = rep(NA, max_iter) ; s_dual10value = rep(NA, max_iter)
  r_primal11value = rep(NA, max_iter) ; s_dual11value = rep(NA, max_iter)
  r_primal12value = rep(NA, max_iter) ; s_dual12value = rep(NA, max_iter)
  r_primal13value = rep(NA, max_iter) ; s_dual13value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  B2 <- matrix(0, p, 1)
  B3 <- matrix(0, p, 1)
  B4 <- matrix(0, p, 1)
  B5 <- matrix(0, p, 1)
  B6 <- matrix(0, p, 1)
  B7 <- matrix(0, p, 1)
  B8 <- matrix(0, p, 1)
  B9 <- matrix(0, p, 1)
  B10 <- matrix(0, p, 1)
  B11 <- matrix(0, p, 1)
  B12 <- matrix(0, p, 1)
  B13 <- matrix(0, p, 1)
  
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma02 <- matrix(0, p, 1)
  Gamma03 <- matrix(0, p, 1)
  Gamma04 <- matrix(0, p, 1)
  Gamma05 <- matrix(0, p, 1)
  Gamma06 <- matrix(0, p, 1)
  Gamma07 <- matrix(0, p, 1)
  Gamma08 <- matrix(0, p, 1)
  Gamma09 <- matrix(0, p, 1)
  Gamma010 <- matrix(0, p, 1)
  Gamma011 <- matrix(0, p, 1)
  Gamma012 <- matrix(0, p, 1)
  Gamma013 <- matrix(0, p, 1)
  
  Gamma1 <- matrix(0, p, 1)
  Gamma2 <- matrix(0, p, 1)
  Gamma3 <- matrix(0, p, 1)
  Gamma4 <- matrix(0, p, 1)
  Gamma5 <- matrix(0, p, 1)
  Gamma6 <- matrix(0, p, 1)
  Gamma7 <- matrix(0, p, 1)
  Gamma8 <- matrix(0, p, 1)
  Gamma9 <- matrix(0, p, 1)
  Gamma10 <- matrix(0, p, 1)
  Gamma11 <- matrix(0, p, 1)
  Gamma12 <- matrix(0, p, 1)
  Gamma13 <- matrix(0, p, 1)
  
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  U2 <- matrix(0, p, 1)
  U3 <- matrix(0, p, 1)
  U4 <- matrix(0, p, 1)
  U5 <- matrix(0, p, 1)
  U6 <- matrix(0, p, 1)
  U7 <- matrix(0, p, 1)
  U8 <- matrix(0, p, 1)
  U9 <- matrix(0, p, 1)
  U10 <- matrix(0, p, 1)
  U11 <- matrix(0, p, 1)
  U12 <- matrix(0, p, 1)
  U13 <- matrix(0, p, 1)
  
  V1 <- matrix(0, p, 1)
  V2 <- matrix(0, p, 1)
  V3 <- matrix(0, p, 1)
  V4 <- matrix(0, p, 1)
  V5 <- matrix(0, p, 1)
  V6 <- matrix(0, p, 1)
  V7 <- matrix(0, p, 1)
  V8 <- matrix(0, p, 1)
  V9 <- matrix(0, p, 1)
  V10 <- matrix(0, p, 1)
  V11 <- matrix(0, p, 1)
  V12 <- matrix(0, p, 1)
  V13 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0)
    X1tX1 = crossprod(X1,X1)
    X2tX2 = crossprod(X2,X2)
    X3tX3 = crossprod(X3,X3)
    X4tX4 = crossprod(X4,X4)
    X5tX5 = crossprod(X5,X5)
    X6tX6 = crossprod(X6,X6)
    X7tX7 = crossprod(X7,X7)
    X8tX8 = crossprod(X8,X8)
    X9tX9 = crossprod(X9,X9)
    X10tX10 = crossprod(X10,X10)
    X11tX11 = crossprod(X11,X11)
    X12tX12 = crossprod(X12,X12)
    X13tX13 = crossprod(X13,X13)
    X0ty0 = crossprod(X0,y0)
    X1ty1 = crossprod(X1,y1)
    X2ty2 = crossprod(X2,y2)
    X3ty3 = crossprod(X3,y3)
    X4ty4 = crossprod(X4,y4)
    X5ty5 = crossprod(X5,y5)
    X6ty6 = crossprod(X6,y6)
    X7ty7 = crossprod(X7,y7)
    X8ty8 = crossprod(X8,y8)
    X9ty9 = crossprod(X9,y9)
    X10ty10 = crossprod(X10,y10)
    X11ty11 = crossprod(X11,y11)
    X12ty12 = crossprod(X12,y12)
    X13ty13 = crossprod(X13,y13)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma03_old = Gamma03
    Gamma04_old = Gamma04
    Gamma05_old = Gamma05
    Gamma06_old = Gamma06
    Gamma07_old = Gamma07
    Gamma08_old = Gamma08
    Gamma09_old = Gamma09
    Gamma010_old = Gamma010
    Gamma011_old = Gamma011
    Gamma012_old = Gamma012
    Gamma013_old = Gamma013
    
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    Gamma3_old = Gamma3
    Gamma4_old = Gamma4
    Gamma5_old = Gamma5
    Gamma6_old = Gamma6
    Gamma7_old = Gamma7
    Gamma8_old = Gamma8
    Gamma9_old = Gamma9
    Gamma10_old = Gamma10
    Gamma11_old = Gamma11
    Gamma12_old = Gamma12
    Gamma13_old = Gamma13
    
    #### Step 1: B0 update ####
    lhs = 1/N*X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06 + rho07 + rho08 + rho09 + rho010 + rho011 + rho012 + rho013)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6) + rho07 * (Gamma07 - U7) + rho08 * (Gamma08 - U8) + rho09 * (Gamma09 - U9) + rho010 * (Gamma010 - U10) + rho011 * (Gamma011 - U11) + rho012 * (Gamma012 - U12) + rho013 * (Gamma013 - U13)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/N*X2tX2 + rho2*Ip
    rhs = 1/N*X2ty2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/N*X3tX3 + rho3*Ip
    rhs = 1/N*X3ty3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/N*X4tX4 + rho4*Ip
    rhs = 1/N*X4ty4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### B5 update ####
    lhs = 1/N*X5tX5 + rho5*Ip
    rhs = 1/N*X5ty5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    #### B6 update ####
    lhs = 1/N*X6tX6 + rho6*Ip
    rhs = 1/N*X6ty6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    
    #### B7 update ####
    lhs = 1/N*X7tX7 + rho7*Ip
    rhs = 1/N*X7ty7 + rho7*(Gamma7 - V7)
    B7 <- solve(lhs, rhs)
    
    #### B8 update ####
    lhs = 1/N*X8tX8 + rho8*Ip
    rhs = 1/N*X8ty8 + rho8*(Gamma8 - V8)
    B8 <- solve(lhs, rhs)
    
    #### B9 update ####
    lhs = 1/N*X9tX9 + rho9*Ip
    rhs = 1/N*X9ty9 + rho9*(Gamma9 - V9)
    B9 <- solve(lhs, rhs)
    
    #### B10 update ####
    lhs = 1/N*X10tX10 + rho10*Ip
    rhs = 1/N*X10ty10 + rho10*(Gamma10 - V10)
    B10 <- solve(lhs, rhs)
    
    #### B11 update ####
    lhs = 1/N*X11tX11 + rho11*Ip
    rhs = 1/N*X11ty11 + rho11*(Gamma11 - V11)
    B11 <- solve(lhs, rhs)
    
    #### B12 update ####
    lhs = 1/N*X12tX12 + rho12*Ip
    rhs = 1/N*X12ty12 + rho12*(Gamma12 - V12)
    B12 <- solve(lhs, rhs)
    
    #### B13 update ####
    lhs = 1/N*X13tX13 + rho13*Ip
    rhs = 1/N*X13ty13 + rho13*(Gamma13 - V13)
    B13 <- solve(lhs, rhs)
    
    #### Gamma00 update ####
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Gamma01 update ####
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    #### Gamma02 update ####
    Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
    #### Gamma03 update ####
    Gamma03 = Gamma3 + soft_threshold(B0+U3-Gamma3, lambda3/rho03)
    #### Gamma04 update ####
    Gamma04 = Gamma4 + soft_threshold(B0+U4-Gamma4, lambda4/rho04)
    #### Gamma05 update ####
    Gamma05 = Gamma5 + soft_threshold(B0+U5-Gamma5, lambda5/rho05)
    #### Gamma06 update ####
    Gamma06 = Gamma6 + soft_threshold(B0+U6-Gamma6, lambda6/rho06)
    #### Gamma07 update ####
    Gamma07 = Gamma7 + soft_threshold(B0+U7-Gamma7, lambda7/rho07)
    #### Gamma08 update ####
    Gamma08 = Gamma8 + soft_threshold(B0+U8-Gamma8, lambda8/rho08)
    #### Gamma09 update ####
    Gamma09 = Gamma9 + soft_threshold(B0+U9-Gamma9, lambda9/rho09)
    #### Gamma010 update ####
    Gamma010 = Gamma10 + soft_threshold(B0+U10-Gamma10, lambda10/rho010)
    #### Gamma011 update ####
    Gamma011 = Gamma11 + soft_threshold(B0+U11-Gamma11, lambda11/rho011)
    #### Gamma012 update ####
    Gamma012 = Gamma12 + soft_threshold(B0+U12-Gamma12, lambda12/rho012)
    #### Gamma013 update ####
    Gamma013 = Gamma13 + soft_threshold(B0+U13-Gamma13, lambda13/rho013)
    
    #### Gamma1 update ####
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
    #### Gamma2 update ####
    Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
    #### Gamma3 update ####
    Gamma3 = Gamma03 + soft_threshold(B3+V3-Gamma03, lambda3/rho3)
    #### Gamma4 update ####
    Gamma4 = Gamma04 + soft_threshold(B4+V4-Gamma04, lambda4/rho4)
    #### Gamma5 update ####
    Gamma5 = Gamma05 + soft_threshold(B5+V5-Gamma05, lambda5/rho5)
    #### Gamma6 update ####
    Gamma6 = Gamma06 + soft_threshold(B6+V6-Gamma06, lambda6/rho6)
    #### Gamma7 update ####
    Gamma7 = Gamma07 + soft_threshold(B7+V7-Gamma07, lambda7/rho7)
    #### Gamma8 update ####
    Gamma8 = Gamma08 + soft_threshold(B8+V8-Gamma08, lambda8/rho8)
    #### Gamma9 update ####
    Gamma9 = Gamma09 + soft_threshold(B9+V9-Gamma09, lambda9/rho9)
    #### Gamma10 update ####
    Gamma10 = Gamma010 + soft_threshold(B10+V10-Gamma010, lambda10/rho10)
    #### Gamma11 update ####
    Gamma11 = Gamma011 + soft_threshold(B11+V11-Gamma011, lambda11/rho11)
    #### Gamma12 update ####
    Gamma12 = Gamma012 + soft_threshold(B12+V12-Gamma012, lambda12/rho12)
    #### Gamma13 update ####
    Gamma13 = Gamma013 + soft_threshold(B13+V13-Gamma013, lambda13/rho13)
    
    #### Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    U3 <- U3 + B0 - Gamma03
    U4 <- U4 + B0 - Gamma04
    U5 <- U5 + B0 - Gamma05
    U6 <- U6 + B0 - Gamma06
    U7 <- U7 + B0 - Gamma07
    U8 <- U8 + B0 - Gamma08
    U9 <- U9 + B0 - Gamma09
    U10 <- U10 + B0 - Gamma010
    U11 <- U11 + B0 - Gamma011
    U12 <- U12 + B0 - Gamma012
    U13 <- U13 + B0 - Gamma013
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    V3 <- V3 + B3 - Gamma3
    V4 <- V4 + B4 - Gamma4
    V5 <- V5 + B5 - Gamma5
    V6 <- V6 + B6 - Gamma6
    V7 <- V7 + B7 - Gamma7
    V8 <- V8 + B8 - Gamma8
    V9 <- V9 + B9 - Gamma9
    V10 <- V10 + B10 - Gamma10
    V11 <- V11 + B11 - Gamma11
    V12 <- V12 + B12 - Gamma12
    V13 <- V13 + B13 - Gamma13
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
    r_primal03 = l2norm(B0 - Gamma03) ; s_dual03 = rho03 * l2norm(Gamma03 - Gamma03_old)
    r_primal04 = l2norm(B0 - Gamma04) ; s_dual04 = rho04 * l2norm(Gamma04 - Gamma04_old)
    r_primal05 = l2norm(B0 - Gamma05) ; s_dual05 = rho05 * l2norm(Gamma05 - Gamma05_old)
    r_primal06 = l2norm(B0 - Gamma06) ; s_dual06 = rho06 * l2norm(Gamma06 - Gamma06_old)
    r_primal07 = l2norm(B0 - Gamma07) ; s_dual07 = rho07 * l2norm(Gamma07 - Gamma07_old)
    r_primal08 = l2norm(B0 - Gamma08) ; s_dual08 = rho08 * l2norm(Gamma08 - Gamma08_old)
    r_primal09 = l2norm(B0 - Gamma09) ; s_dual09 = rho09 * l2norm(Gamma09 - Gamma09_old)
    r_primal010 = l2norm(B0 - Gamma010) ; s_dual010 = rho010 * l2norm(Gamma010 - Gamma010_old)
    r_primal011 = l2norm(B0 - Gamma011) ; s_dual011 = rho011 * l2norm(Gamma011 - Gamma011_old)
    r_primal012 = l2norm(B0 - Gamma012) ; s_dual012 = rho012 * l2norm(Gamma012 - Gamma012_old)
    r_primal013 = l2norm(B0 - Gamma013) ; s_dual013 = rho013 * l2norm(Gamma013 - Gamma013_old)
    
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
    r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
    r_primal3 = l2norm(B3 - Gamma3) ; s_dual3 = rho3 * l2norm(Gamma3 - Gamma3_old)
    r_primal4 = l2norm(B4 - Gamma4) ; s_dual4 = rho4 * l2norm(Gamma4 - Gamma4_old)
    r_primal5 = l2norm(B5 - Gamma5) ; s_dual5 = rho5 * l2norm(Gamma5 - Gamma5_old)
    r_primal6 = l2norm(B6 - Gamma6) ; s_dual6 = rho6 * l2norm(Gamma6 - Gamma6_old)
    r_primal7 = l2norm(B7 - Gamma7) ; s_dual7 = rho7 * l2norm(Gamma7 - Gamma7_old)
    r_primal8 = l2norm(B8 - Gamma8) ; s_dual8 = rho8 * l2norm(Gamma8 - Gamma8_old)
    r_primal9 = l2norm(B9 - Gamma9) ; s_dual9 = rho9 * l2norm(Gamma9 - Gamma9_old)
    r_primal10 = l2norm(B10 - Gamma10) ; s_dual10 = rho10 * l2norm(Gamma10 - Gamma10_old)
    r_primal11 = l2norm(B11 - Gamma11) ; s_dual11 = rho11 * l2norm(Gamma11 - Gamma11_old)
    r_primal12 = l2norm(B12 - Gamma12) ; s_dual12 = rho12 * l2norm(Gamma12 - Gamma12_old)
    r_primal13 = l2norm(B13 - Gamma13) ; s_dual13 = rho13 * l2norm(Gamma13 - Gamma13_old)
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
    r_primal03value[iter] = r_primal03 ; s_dual03value[iter] = s_dual03
    r_primal04value[iter] = r_primal04 ; s_dual04value[iter] = s_dual04
    r_primal05value[iter] = r_primal05 ; s_dual05value[iter] = s_dual05
    r_primal06value[iter] = r_primal06 ; s_dual06value[iter] = s_dual06
    r_primal07value[iter] = r_primal07 ; s_dual07value[iter] = s_dual07
    r_primal08value[iter] = r_primal08 ; s_dual08value[iter] = s_dual08
    r_primal09value[iter] = r_primal09 ; s_dual09value[iter] = s_dual09
    r_primal010value[iter] = r_primal010 ; s_dual010value[iter] = s_dual010
    r_primal011value[iter] = r_primal011 ; s_dual011value[iter] = s_dual011
    r_primal012value[iter] = r_primal012 ; s_dual012value[iter] = s_dual012
    r_primal013value[iter] = r_primal013 ; s_dual013value[iter] = s_dual013
    
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
    r_primal3value[iter] = r_primal3 ; s_dual3value[iter] = s_dual3
    r_primal4value[iter] = r_primal4 ; s_dual4value[iter] = s_dual4
    r_primal5value[iter] = r_primal5 ; s_dual5value[iter] = s_dual5
    r_primal6value[iter] = r_primal6 ; s_dual6value[iter] = s_dual6
    r_primal7value[iter] = r_primal7 ; s_dual7value[iter] = s_dual7
    r_primal8value[iter] = r_primal8 ; s_dual8value[iter] = s_dual8
    r_primal9value[iter] = r_primal9 ; s_dual9value[iter] = s_dual9
    r_primal10value[iter] = r_primal10 ; s_dual10value[iter] = s_dual10
    r_primal11value[iter] = r_primal11 ; s_dual11value[iter] = s_dual11
    r_primal12value[iter] = r_primal12 ; s_dual12value[iter] = s_dual12
    r_primal13value[iter] = r_primal13 ; s_dual13value[iter] = s_dual13
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal02 > mu * s_dual02){
      rho02 = min(rho02 * tau_incr, rho_max)
      U2 = U2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual02 > mu * r_primal02){
      rho02 = rho02 / tau_decr
      U2 = U2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal03 > mu * s_dual03){
      rho03 = min(rho03 * tau_incr, rho_max)
      U3 = U3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual03 > mu * r_primal03){
      rho03 = rho03 / tau_decr
      U3 = U3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal04 > mu * s_dual04){
      rho04 = min(rho04 * tau_incr, rho_max)
      U4 = U4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual04 > mu * r_primal04){
      rho04 = rho04 / tau_decr
      U4 = U4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal05 > mu * s_dual05){
      rho05 = min(rho05 * tau_incr, rho_max)
      U5 = U5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual05 > mu * r_primal05){
      rho05 = rho05 / tau_decr
      U5 = U5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal06 > mu * s_dual06){
      rho06 = min(rho06 * tau_incr, rho_max)
      U6 = U6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual06 > mu * r_primal06){
      rho06 = rho06 / tau_decr
      U6 = U6 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal07 > mu * s_dual07){
      rho07 = min(rho07 * tau_incr, rho_max)
      U7 = U7 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual07 > mu * r_primal07){
      rho07 = rho07 / tau_decr
      U7 = U7 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal08 > mu * s_dual08){
      rho08 = min(rho08 * tau_incr, rho_max)
      U8 = U8 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual08 > mu * r_primal08){
      rho08 = rho08 / tau_decr
      U8 = U8 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal09 > mu * s_dual09){
      rho09 = min(rho09 * tau_incr, rho_max)
      U9 = U9 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual09 > mu * r_primal09){
      rho09 = rho09 / tau_decr
      U9 = U9 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal010 > mu * s_dual010){
      rho010 = min(rho010 * tau_incr, rho_max)
      U10 = U10 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual010 > mu * r_primal010){
      rho010 = rho010 / tau_decr
      U10 = U10 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal011 > mu * s_dual011){
      rho011 = min(rho011 * tau_incr, rho_max)
      U11 = U11 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual011 > mu * r_primal011){
      rho011 = rho011 / tau_decr
      U11 = U11 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal012 > mu * s_dual012){
      rho012 = min(rho012 * tau_incr, rho_max)
      U12 = U12 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual012 > mu * r_primal012){
      rho012 = rho012 / tau_decr
      U12 = U12 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal013 > mu * s_dual013){
      rho013 = min(rho013 * tau_incr, rho_max)
      U13 = U13 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual013 > mu * r_primal013){
      rho013 = rho013 / tau_decr
      U13 = U13 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal2 > mu * s_dual2){
      rho2 = min(rho2 * tau_incr, rho_max)
      V2 = V2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual2 > mu * r_primal2){
      rho2 = rho2 / tau_decr
      V2 = V2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal3 > mu * s_dual3){
      rho3 = min(rho3 * tau_incr, rho_max)
      V3 = V3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual3 > mu * r_primal3){
      rho3 = rho3 / tau_decr
      V3 = V3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal4 > mu * s_dual4){
      rho4 = min(rho4 * tau_incr, rho_max)
      V4 = V4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual4 > mu * r_primal4){
      rho4 = rho4 / tau_decr
      V4 = V4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal5 > mu * s_dual5){
      rho5 = min(rho5 * tau_incr, rho_max)
      V5 = V5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual5 > mu * r_primal5){
      rho5 = rho5 / tau_decr
      V5 = V5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal6 > mu * s_dual6){
      rho6 = min(rho6 * tau_incr, rho_max)
      V6 = V6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual6 > mu * r_primal6){
      rho6 = rho6 / tau_decr
      V6 = V6 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal7 > mu * s_dual7){
      rho7 = min(rho7 * tau_incr, rho_max)
      V7 = V7 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual7 > mu * r_primal7){
      rho7 = rho7 / tau_decr
      V7 = V7 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal8 > mu * s_dual8){
      rho8 = min(rho8 * tau_incr, rho_max)
      V8 = V8 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual8 > mu * r_primal8){
      rho8 = rho8 / tau_decr
      V8 = V8 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal9 > mu * s_dual9){
      rho9 = min(rho9 * tau_incr, rho_max)
      V9 = V9 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual9 > mu * r_primal9){
      rho9 = rho9 / tau_decr
      V9 = V9 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal10 > mu * s_dual10){
      rho10 = min(rho10 * tau_incr, rho_max)
      V10 = V10 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual10 > mu * r_primal10){
      rho10 = rho10 / tau_decr
      V10 = V10 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal11 > mu * s_dual11){
      rho11 = min(rho11 * tau_incr, rho_max)
      V11 = V11 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual11 > mu * r_primal11){
      rho11 = rho11 / tau_decr
      V11 = V11 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal12 > mu * s_dual12){
      rho12 = min(rho12 * tau_incr, rho_max)
      V12 = V12 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual12 > mu * r_primal12){
      rho12 = rho12 / tau_decr
      V12 = V12 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal13 > mu * s_dual13){
      rho13 = min(rho13 * tau_incr, rho_max)
      V13 = V13 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual13 > mu * r_primal13){
      rho13 = rho13 / tau_decr
      V13 = V13 * tau_decr # change of variable required for the scaled variable
    }
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho02value[iter] = rho02
    rho03value[iter] = rho03
    rho04value[iter] = rho04
    rho05value[iter] = rho05
    rho06value[iter] = rho06
    rho07value[iter] = rho07
    rho08value[iter] = rho08
    rho09value[iter] = rho09
    rho010value[iter] = rho010
    rho011value[iter] = rho011
    rho012value[iter] = rho012
    rho013value[iter] = rho013
    
    rho1value[iter] = rho1
    rho2value[iter] = rho2
    rho3value[iter] = rho3
    rho4value[iter] = rho4
    rho5value[iter] = rho5
    rho6value[iter] = rho6
    rho7value[iter] = rho7
    rho8value[iter] = rho8
    rho9value[iter] = rho9
    rho10value[iter] = rho10
    rho11value[iter] = rho11
    rho12value[iter] = rho12
    rho13value[iter] = rho13
    
    if (max(r_primal00,
            r_primal01, r_primal02, r_primal03, r_primal04, r_primal05, r_primal06, r_primal07, r_primal08, r_primal09, r_primal010, r_primal011, r_primal012, r_primal013,
            r_primal1, r_primal2, r_primal3, r_primal4, r_primal5, r_primal6, r_primal7, r_primal8, r_primal9, r_primal10, r_primal11, r_primal12, r_primal13) < tol_prim
        &&
        max(s_dual00,
            s_dual01, s_dual02, s_dual03, s_dual04, s_dual05, s_dual06, s_dual07, s_dual08, s_dual09, s_dual010, s_dual011, s_dual012, s_dual013,
            s_dual1, s_dual2, s_dual3, s_dual4, s_dual5, s_dual6, s_dual7, s_dual8, s_dual9, s_dual10, s_dual11, s_dual12, s_dual13) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6, B7 = B7, B8 = B8, B9 = B9, B10 = B10, B11 = B11, B12 = B12, B13 = B13,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06, Gamma07 = Gamma07, Gamma08 = Gamma08, Gamma09 = Gamma09, Gamma010 = Gamma010, Gamma011 = Gamma011, Gamma012 = Gamma012, Gamma013 = Gamma013,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6, Gamma7 = Gamma7, Gamma8 = Gamma8, Gamma9 = Gamma9, Gamma10 = Gamma10, Gamma11 = Gamma11, Gamma12 = Gamma12, Gamma13 = Gamma13,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6, U7 = U7, U8 = U8, U9 = U9, U10 = U10, U11 = U11, U12 = U12, U13 = U13,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6, V7 = V7, V8 = V8, V9 = V9, V10 = V10, V11 = V11, V12 = V12, V13 = V13,
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, rho04 = rho04value, rho05 = rho05value, rho06 = rho06value, rho07 = rho07value, rho08 = rho08value, rho09 = rho09value, rho010 = rho010value, rho011 = rho011value, rho012 = rho012value, rho013 = rho013value,
              rho1 = rho1value, rho2 = rho2value, rho3 = rho3value, rho4 = rho4value, rho5 = rho5value, rho6 = rho6value, rho7 = rho7value, rho8 = rho8value, rho9 = rho9value, rho10 = rho10value, rho11 = rho11value, rho12 = rho12value, rho13 = rho13value,
              
              r_primal00 = r_primal00value, s_dual00 = s_dual00value,
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              r_primal04 = r_primal04value, s_dual04 = s_dual04value,
              r_primal05 = r_primal05value, s_dual05 = s_dual05value,
              r_primal06 = r_primal06value, s_dual06 = s_dual06value,
              r_primal07 = r_primal07value, s_dual07 = s_dual07value,
              r_primal08 = r_primal08value, s_dual08 = s_dual08value,
              r_primal09 = r_primal09value, s_dual09 = s_dual09value,
              r_primal010 = r_primal010value, s_dual010 = s_dual010value,
              r_primal011 = r_primal011value, s_dual011 = s_dual011value,
              r_primal012 = r_primal012value, s_dual012 = s_dual012value,
              r_primal013 = r_primal013value, s_dual013 = s_dual013value,
              
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value,
              r_primal4 = r_primal4value, s_dual4 = s_dual4value,
              r_primal5 = r_primal5value, s_dual5 = s_dual5value,
              r_primal6 = r_primal6value, s_dual6 = s_dual6value,
              r_primal7 = r_primal7value, s_dual7 = s_dual7value,
              r_primal8 = r_primal8value, s_dual8 = s_dual8value,
              r_primal9 = r_primal9value, s_dual9 = s_dual9value,
              r_primal10 = r_primal10value, s_dual10 = s_dual10value,
              r_primal11 = r_primal11value, s_dual11 = s_dual11value,
              r_primal12 = r_primal12value, s_dual12 = s_dual12value,
              r_primal13 = r_primal13value, s_dual13 = s_dual13value))
}

admm_stl_ninesource <- function(X0, y0, X1, y1, X2, y2, X3, y3, X4, y4, X5, y5, X6, y6, X7, y7, X8, y8, X9, y9,
                                X0tX0 = NULL, X0ty0 = NULL,
                                X1tX1 = NULL, X1ty1 = NULL,
                                X2tX2 = NULL, X2ty2 = NULL,
                                X3tX3 = NULL, X3ty3 = NULL,
                                X4tX4 = NULL, X4ty4 = NULL,
                                X5tX5 = NULL, X5ty5 = NULL,
                                X6tX6 = NULL, X6ty6 = NULL,
                                X7tX7 = NULL, X7ty7 = NULL,
                                X8tX8 = NULL, X8ty8 = NULL,
                                X9tX9 = NULL, X9ty9 = NULL,
                                lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8, lambda9,
                                rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1, rho07 = 1, rho08 = 1, rho09 = 1,
                                rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1, rho7 = 1, rho8 = 1, rho9 = 1,
                                mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                max_iter = 5000, tol_prim = 1e-6, tol_dual = 1e-6, verbose = F){
  
  l2norm = function(x){sqrt(sum(x^2))}
  l2norm_sq = function(x){sum(x^2)}
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3) + nrow(X4) + nrow(X5) + nrow(X6) + nrow(X7) + nrow(X8) + nrow(X9)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho02value = rep(NA, max_iter)
  rho03value = rep(NA, max_iter)
  rho04value = rep(NA, max_iter)
  rho05value = rep(NA, max_iter)
  rho06value = rep(NA, max_iter)
  rho07value = rep(NA, max_iter)
  rho08value = rep(NA, max_iter)
  rho09value = rep(NA, max_iter)
  
  rho1value = rep(NA, max_iter)
  rho2value = rep(NA, max_iter)
  rho3value = rep(NA, max_iter)
  rho4value = rep(NA, max_iter)
  rho5value = rep(NA, max_iter)
  rho6value = rep(NA, max_iter)
  rho7value = rep(NA, max_iter)
  rho8value = rep(NA, max_iter)
  rho9value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
  r_primal04value = rep(NA, max_iter) ; s_dual04value = rep(NA, max_iter)
  r_primal05value = rep(NA, max_iter) ; s_dual05value = rep(NA, max_iter)
  r_primal06value = rep(NA, max_iter) ; s_dual06value = rep(NA, max_iter)
  r_primal07value = rep(NA, max_iter) ; s_dual07value = rep(NA, max_iter)
  r_primal08value = rep(NA, max_iter) ; s_dual08value = rep(NA, max_iter)
  r_primal09value = rep(NA, max_iter) ; s_dual09value = rep(NA, max_iter)
  
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
  r_primal4value = rep(NA, max_iter) ; s_dual4value = rep(NA, max_iter)
  r_primal5value = rep(NA, max_iter) ; s_dual5value = rep(NA, max_iter)
  r_primal6value = rep(NA, max_iter) ; s_dual6value = rep(NA, max_iter)
  r_primal7value = rep(NA, max_iter) ; s_dual7value = rep(NA, max_iter)
  r_primal8value = rep(NA, max_iter) ; s_dual8value = rep(NA, max_iter)
  r_primal9value = rep(NA, max_iter) ; s_dual9value = rep(NA, max_iter)
  
  p <- ncol(X0)
  
  B0 <- matrix(0, p, 1)
  B1 <- matrix(0, p, 1)
  B2 <- matrix(0, p, 1)
  B3 <- matrix(0, p, 1)
  B4 <- matrix(0, p, 1)
  B5 <- matrix(0, p, 1)
  B6 <- matrix(0, p, 1)
  B7 <- matrix(0, p, 1)
  B8 <- matrix(0, p, 1)
  B9 <- matrix(0, p, 1)
  
  Gamma00 <- matrix(0, p, 1)
  Gamma01 <- matrix(0, p, 1)
  Gamma02 <- matrix(0, p, 1)
  Gamma03 <- matrix(0, p, 1)
  Gamma04 <- matrix(0, p, 1)
  Gamma05 <- matrix(0, p, 1)
  Gamma06 <- matrix(0, p, 1)
  Gamma07 <- matrix(0, p, 1)
  Gamma08 <- matrix(0, p, 1)
  Gamma09 <- matrix(0, p, 1)
  
  Gamma1 <- matrix(0, p, 1)
  Gamma2 <- matrix(0, p, 1)
  Gamma3 <- matrix(0, p, 1)
  Gamma4 <- matrix(0, p, 1)
  Gamma5 <- matrix(0, p, 1)
  Gamma6 <- matrix(0, p, 1)
  Gamma7 <- matrix(0, p, 1)
  Gamma8 <- matrix(0, p, 1)
  Gamma9 <- matrix(0, p, 1)
  
  U0 <- matrix(0, p, 1)
  U1 <- matrix(0, p, 1)
  U2 <- matrix(0, p, 1)
  U3 <- matrix(0, p, 1)
  U4 <- matrix(0, p, 1)
  U5 <- matrix(0, p, 1)
  U6 <- matrix(0, p, 1)
  U7 <- matrix(0, p, 1)
  U8 <- matrix(0, p, 1)
  U9 <- matrix(0, p, 1)
  
  V1 <- matrix(0, p, 1)
  V2 <- matrix(0, p, 1)
  V3 <- matrix(0, p, 1)
  V4 <- matrix(0, p, 1)
  V5 <- matrix(0, p, 1)
  V6 <- matrix(0, p, 1)
  V7 <- matrix(0, p, 1)
  V8 <- matrix(0, p, 1)
  V9 <- matrix(0, p, 1)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0)
    X1tX1 = crossprod(X1,X1)
    X2tX2 = crossprod(X2,X2)
    X3tX3 = crossprod(X3,X3)
    X4tX4 = crossprod(X4,X4)
    X5tX5 = crossprod(X5,X5)
    X6tX6 = crossprod(X6,X6)
    X7tX7 = crossprod(X7,X7)
    X8tX8 = crossprod(X8,X8)
    X9tX9 = crossprod(X9,X9)
    X0ty0 = crossprod(X0,y0)
    X1ty1 = crossprod(X1,y1)
    X2ty2 = crossprod(X2,y2)
    X3ty3 = crossprod(X3,y3)
    X4ty4 = crossprod(X4,y4)
    X5ty5 = crossprod(X5,y5)
    X6ty6 = crossprod(X6,y6)
    X7ty7 = crossprod(X7,y7)
    X8ty8 = crossprod(X8,y8)
    X9ty9 = crossprod(X9,y9)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma03_old = Gamma03
    Gamma04_old = Gamma04
    Gamma05_old = Gamma05
    Gamma06_old = Gamma06
    Gamma07_old = Gamma07
    Gamma08_old = Gamma08
    Gamma09_old = Gamma09
    
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    Gamma3_old = Gamma3
    Gamma4_old = Gamma4
    Gamma5_old = Gamma5
    Gamma6_old = Gamma6
    Gamma7_old = Gamma7
    Gamma8_old = Gamma8
    Gamma9_old = Gamma9
    
    #### Step 1: B0 update ####
    lhs = 1/N*X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06 + rho07 + rho08 + rho09)*Ip
    rhs = 1/N*X0ty0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6) + rho07 * (Gamma07 - U7) + rho08 * (Gamma08 - U8) + rho09 * (Gamma09 - U9)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/N*X1tX1 + rho1*Ip
    rhs = 1/N*X1ty1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/N*X2tX2 + rho2*Ip
    rhs = 1/N*X2ty2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/N*X3tX3 + rho3*Ip
    rhs = 1/N*X3ty3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/N*X4tX4 + rho4*Ip
    rhs = 1/N*X4ty4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### B5 update ####
    lhs = 1/N*X5tX5 + rho5*Ip
    rhs = 1/N*X5ty5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    #### B6 update ####
    lhs = 1/N*X6tX6 + rho6*Ip
    rhs = 1/N*X6ty6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    
    #### B7 update ####
    lhs = 1/N*X7tX7 + rho7*Ip
    rhs = 1/N*X7ty7 + rho7*(Gamma7 - V7)
    B7 <- solve(lhs, rhs)
    
    #### B8 update ####
    lhs = 1/N*X8tX8 + rho8*Ip
    rhs = 1/N*X8ty8 + rho8*(Gamma8 - V8)
    B8 <- solve(lhs, rhs)
    
    #### B9 update ####
    lhs = 1/N*X9tX9 + rho9*Ip
    rhs = 1/N*X9ty9 + rho9*(Gamma9 - V9)
    B9 <- solve(lhs, rhs)
    
    #### Gamma00 update ####
    Gamma00 = soft_threshold(B0+U0, lambda0/rho00)
    
    #### Gamma01 update ####
    Gamma01 = Gamma1 + soft_threshold(B0+U1-Gamma1, lambda1/rho01)
    #### Gamma02 update ####
    Gamma02 = Gamma2 + soft_threshold(B0+U2-Gamma2, lambda2/rho02)
    #### Gamma03 update ####
    Gamma03 = Gamma3 + soft_threshold(B0+U3-Gamma3, lambda3/rho03)
    #### Gamma04 update ####
    Gamma04 = Gamma4 + soft_threshold(B0+U4-Gamma4, lambda4/rho04)
    #### Gamma05 update ####
    Gamma05 = Gamma5 + soft_threshold(B0+U5-Gamma5, lambda5/rho05)
    #### Gamma06 update ####
    Gamma06 = Gamma6 + soft_threshold(B0+U6-Gamma6, lambda6/rho06)
    #### Gamma07 update ####
    Gamma07 = Gamma7 + soft_threshold(B0+U7-Gamma7, lambda7/rho07)
    #### Gamma08 update ####
    Gamma08 = Gamma8 + soft_threshold(B0+U8-Gamma8, lambda8/rho08)
    #### Gamma09 update ####
    Gamma09 = Gamma9 + soft_threshold(B0+U9-Gamma9, lambda9/rho09)
    
    #### Gamma1 update ####
    Gamma1 = Gamma01 + soft_threshold(B1+V1-Gamma01, lambda1/rho1)
    #### Gamma2 update ####
    Gamma2 = Gamma02 + soft_threshold(B2+V2-Gamma02, lambda2/rho2)
    #### Gamma3 update ####
    Gamma3 = Gamma03 + soft_threshold(B3+V3-Gamma03, lambda3/rho3)
    #### Gamma4 update ####
    Gamma4 = Gamma04 + soft_threshold(B4+V4-Gamma04, lambda4/rho4)
    #### Gamma5 update ####
    Gamma5 = Gamma05 + soft_threshold(B5+V5-Gamma05, lambda5/rho5)
    #### Gamma6 update ####
    Gamma6 = Gamma06 + soft_threshold(B6+V6-Gamma06, lambda6/rho6)
    #### Gamma7 update ####
    Gamma7 = Gamma07 + soft_threshold(B7+V7-Gamma07, lambda7/rho7)
    #### Gamma8 update ####
    Gamma8 = Gamma08 + soft_threshold(B8+V8-Gamma08, lambda8/rho8)
    #### Gamma9 update ####
    Gamma9 = Gamma09 + soft_threshold(B9+V9-Gamma09, lambda9/rho9)
    
    #### Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    U3 <- U3 + B0 - Gamma03
    U4 <- U4 + B0 - Gamma04
    U5 <- U5 + B0 - Gamma05
    U6 <- U6 + B0 - Gamma06
    U7 <- U7 + B0 - Gamma07
    U8 <- U8 + B0 - Gamma08
    U9 <- U9 + B0 - Gamma09
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    V3 <- V3 + B3 - Gamma3
    V4 <- V4 + B4 - Gamma4
    V5 <- V5 + B5 - Gamma5
    V6 <- V6 + B6 - Gamma6
    V7 <- V7 + B7 - Gamma7
    V8 <- V8 + B8 - Gamma8
    V9 <- V9 + B9 - Gamma9
    
    #### Convergence check: primal & dual ####
    r_primal00 = l2norm(B0 - Gamma00) ; s_dual00 = rho00 * l2norm(Gamma00 - Gamma00_old)
    r_primal01 = l2norm(B0 - Gamma01) ; s_dual01 = rho01 * l2norm(Gamma01 - Gamma01_old)
    r_primal02 = l2norm(B0 - Gamma02) ; s_dual02 = rho02 * l2norm(Gamma02 - Gamma02_old)
    r_primal03 = l2norm(B0 - Gamma03) ; s_dual03 = rho03 * l2norm(Gamma03 - Gamma03_old)
    r_primal04 = l2norm(B0 - Gamma04) ; s_dual04 = rho04 * l2norm(Gamma04 - Gamma04_old)
    r_primal05 = l2norm(B0 - Gamma05) ; s_dual05 = rho05 * l2norm(Gamma05 - Gamma05_old)
    r_primal06 = l2norm(B0 - Gamma06) ; s_dual06 = rho06 * l2norm(Gamma06 - Gamma06_old)
    r_primal07 = l2norm(B0 - Gamma07) ; s_dual07 = rho07 * l2norm(Gamma07 - Gamma07_old)
    r_primal08 = l2norm(B0 - Gamma08) ; s_dual08 = rho08 * l2norm(Gamma08 - Gamma08_old)
    r_primal09 = l2norm(B0 - Gamma09) ; s_dual09 = rho09 * l2norm(Gamma09 - Gamma09_old)
    
    r_primal1 = l2norm(B1 - Gamma1) ; s_dual1 = rho1 * l2norm(Gamma1 - Gamma1_old)
    r_primal2 = l2norm(B2 - Gamma2) ; s_dual2 = rho2 * l2norm(Gamma2 - Gamma2_old)
    r_primal3 = l2norm(B3 - Gamma3) ; s_dual3 = rho3 * l2norm(Gamma3 - Gamma3_old)
    r_primal4 = l2norm(B4 - Gamma4) ; s_dual4 = rho4 * l2norm(Gamma4 - Gamma4_old)
    r_primal5 = l2norm(B5 - Gamma5) ; s_dual5 = rho5 * l2norm(Gamma5 - Gamma5_old)
    r_primal6 = l2norm(B6 - Gamma6) ; s_dual6 = rho6 * l2norm(Gamma6 - Gamma6_old)
    r_primal7 = l2norm(B7 - Gamma7) ; s_dual7 = rho7 * l2norm(Gamma7 - Gamma7_old)
    r_primal8 = l2norm(B8 - Gamma8) ; s_dual8 = rho8 * l2norm(Gamma8 - Gamma8_old)
    r_primal9 = l2norm(B9 - Gamma9) ; s_dual9 = rho9 * l2norm(Gamma9 - Gamma9_old)
    
    r_primal00value[iter] = r_primal00 ; s_dual00value[iter] = s_dual00
    r_primal01value[iter] = r_primal01 ; s_dual01value[iter] = s_dual01
    r_primal02value[iter] = r_primal02 ; s_dual02value[iter] = s_dual02
    r_primal03value[iter] = r_primal03 ; s_dual03value[iter] = s_dual03
    r_primal04value[iter] = r_primal04 ; s_dual04value[iter] = s_dual04
    r_primal05value[iter] = r_primal05 ; s_dual05value[iter] = s_dual05
    r_primal06value[iter] = r_primal06 ; s_dual06value[iter] = s_dual06
    r_primal07value[iter] = r_primal07 ; s_dual07value[iter] = s_dual07
    r_primal08value[iter] = r_primal08 ; s_dual08value[iter] = s_dual08
    r_primal09value[iter] = r_primal09 ; s_dual09value[iter] = s_dual09
    
    r_primal1value[iter] = r_primal1 ; s_dual1value[iter] = s_dual1
    r_primal2value[iter] = r_primal2 ; s_dual2value[iter] = s_dual2
    r_primal3value[iter] = r_primal3 ; s_dual3value[iter] = s_dual3
    r_primal4value[iter] = r_primal4 ; s_dual4value[iter] = s_dual4
    r_primal5value[iter] = r_primal5 ; s_dual5value[iter] = s_dual5
    r_primal6value[iter] = r_primal6 ; s_dual6value[iter] = s_dual6
    r_primal7value[iter] = r_primal7 ; s_dual7value[iter] = s_dual7
    r_primal8value[iter] = r_primal8 ; s_dual8value[iter] = s_dual8
    r_primal9value[iter] = r_primal9 ; s_dual9value[iter] = s_dual9
    
    #Adaptive rho adjustment (Boyd 2011)
    if(r_primal00 > mu * s_dual00){
      rho00 = min(rho00 * tau_incr, rho_max)
      U0 = U0 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual00 > mu * r_primal00){
      rho00 = rho00 / tau_decr
      U0 = U0 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal01 > mu * s_dual01){
      rho01 = min(rho01 * tau_incr, rho_max)
      U1 = U1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual01 > mu * r_primal01){
      rho01 = rho01 / tau_decr
      U1 = U1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal02 > mu * s_dual02){
      rho02 = min(rho02 * tau_incr, rho_max)
      U2 = U2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual02 > mu * r_primal02){
      rho02 = rho02 / tau_decr
      U2 = U2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal03 > mu * s_dual03){
      rho03 = min(rho03 * tau_incr, rho_max)
      U3 = U3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual03 > mu * r_primal03){
      rho03 = rho03 / tau_decr
      U3 = U3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal04 > mu * s_dual04){
      rho04 = min(rho04 * tau_incr, rho_max)
      U4 = U4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual04 > mu * r_primal04){
      rho04 = rho04 / tau_decr
      U4 = U4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal05 > mu * s_dual05){
      rho05 = min(rho05 * tau_incr, rho_max)
      U5 = U5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual05 > mu * r_primal05){
      rho05 = rho05 / tau_decr
      U5 = U5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal06 > mu * s_dual06){
      rho06 = min(rho06 * tau_incr, rho_max)
      U6 = U6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual06 > mu * r_primal06){
      rho06 = rho06 / tau_decr
      U6 = U6 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal07 > mu * s_dual07){
      rho07 = min(rho07 * tau_incr, rho_max)
      U7 = U7 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual07 > mu * r_primal07){
      rho07 = rho07 / tau_decr
      U7 = U7 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal08 > mu * s_dual08){
      rho08 = min(rho08 * tau_incr, rho_max)
      U8 = U8 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual08 > mu * r_primal08){
      rho08 = rho08 / tau_decr
      U8 = U8 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal09 > mu * s_dual09){
      rho09 = min(rho09 * tau_incr, rho_max)
      U9 = U9 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual09 > mu * r_primal09){
      rho09 = rho09 / tau_decr
      U9 = U9 * tau_decr # change of variable required for the scaled variable
    }
    
    if(r_primal1 > mu * s_dual1){
      rho1 = min(rho1 * tau_incr, rho_max)
      V1 = V1 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual1 > mu * r_primal1){
      rho1 = rho1 / tau_decr
      V1 = V1 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal2 > mu * s_dual2){
      rho2 = min(rho2 * tau_incr, rho_max)
      V2 = V2 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual2 > mu * r_primal2){
      rho2 = rho2 / tau_decr
      V2 = V2 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal3 > mu * s_dual3){
      rho3 = min(rho3 * tau_incr, rho_max)
      V3 = V3 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual3 > mu * r_primal3){
      rho3 = rho3 / tau_decr
      V3 = V3 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal4 > mu * s_dual4){
      rho4 = min(rho4 * tau_incr, rho_max)
      V4 = V4 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual4 > mu * r_primal4){
      rho4 = rho4 / tau_decr
      V4 = V4 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal5 > mu * s_dual5){
      rho5 = min(rho5 * tau_incr, rho_max)
      V5 = V5 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual5 > mu * r_primal5){
      rho5 = rho5 / tau_decr
      V5 = V5 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal6 > mu * s_dual6){
      rho6 = min(rho6 * tau_incr, rho_max)
      V6 = V6 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual6 > mu * r_primal6){
      rho6 = rho6 / tau_decr
      V6 = V6 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal7 > mu * s_dual7){
      rho7 = min(rho7 * tau_incr, rho_max)
      V7 = V7 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual7 > mu * r_primal7){
      rho7 = rho7 / tau_decr
      V7 = V7 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal8 > mu * s_dual8){
      rho8 = min(rho8 * tau_incr, rho_max)
      V8 = V8 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual8 > mu * r_primal8){
      rho8 = rho8 / tau_decr
      V8 = V8 * tau_decr # change of variable required for the scaled variable
    }
    if(r_primal9 > mu * s_dual9){
      rho9 = min(rho9 * tau_incr, rho_max)
      V9 = V9 / tau_incr # change of variable required for the scaled variable
    }else if(s_dual9 > mu * r_primal9){
      rho9 = rho9 / tau_decr
      V9 = V9 * tau_decr # change of variable required for the scaled variable
    }
    
    rho00value[iter] = rho00
    rho01value[iter] = rho01
    rho02value[iter] = rho02
    rho03value[iter] = rho03
    rho04value[iter] = rho04
    rho05value[iter] = rho05
    rho06value[iter] = rho06
    rho07value[iter] = rho07
    rho08value[iter] = rho08
    rho09value[iter] = rho09
    
    rho1value[iter] = rho1
    rho2value[iter] = rho2
    rho3value[iter] = rho3
    rho4value[iter] = rho4
    rho5value[iter] = rho5
    rho6value[iter] = rho6
    rho7value[iter] = rho7
    rho8value[iter] = rho8
    rho9value[iter] = rho9
    
    if (max(r_primal00,
            r_primal01, r_primal02, r_primal03, r_primal04, r_primal05, r_primal06, r_primal07, r_primal08, r_primal09,
            r_primal1, r_primal2, r_primal3, r_primal4, r_primal5, r_primal6, r_primal7, r_primal8, r_primal9) < tol_prim
        &&
        max(s_dual00,
            s_dual01, s_dual02, s_dual03, s_dual04, s_dual05, s_dual06, s_dual07, s_dual08, s_dual09,
            s_dual1, s_dual2, s_dual3, s_dual4, s_dual5, s_dual6, s_dual7, s_dual8, s_dual9) < tol_dual) {
      if(verbose){
        cat("ADMM converged at iteration", iter, "\n")
      }
      break
    }
  }
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6, B7 = B7, B8 = B8, B9 = B9,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06, Gamma07 = Gamma07, Gamma08 = Gamma08, Gamma09 = Gamma09,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6, Gamma7 = Gamma7, Gamma8 = Gamma8, Gamma9 = Gamma9,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6, U7 = U7, U8 = U8, U9 = U9,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6, V7 = V7, V8 = V8, V9 = V9,
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, rho04 = rho04value, rho05 = rho05value, rho06 = rho06value, rho07 = rho07value, rho08 = rho08value, rho09 = rho09value,
              rho1 = rho1value, rho2 = rho2value, rho3 = rho3value, rho4 = rho4value, rho5 = rho5value, rho6 = rho6value, rho7 = rho7value, rho8 = rho8value, rho9 = rho9value,
              
              r_primal00 = r_primal00value, s_dual00 = s_dual00value,
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              r_primal04 = r_primal04value, s_dual04 = s_dual04value,
              r_primal05 = r_primal05value, s_dual05 = s_dual05value,
              r_primal06 = r_primal06value, s_dual06 = s_dual06value,
              r_primal07 = r_primal07value, s_dual07 = s_dual07value,
              r_primal08 = r_primal08value, s_dual08 = s_dual08value,
              r_primal09 = r_primal09value, s_dual09 = s_dual09value,
              
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value,
              r_primal4 = r_primal4value, s_dual4 = s_dual4value,
              r_primal5 = r_primal5value, s_dual5 = s_dual5value,
              r_primal6 = r_primal6value, s_dual6 = s_dual6value,
              r_primal7 = r_primal7value, s_dual7 = s_dual7value,
              r_primal8 = r_primal8value, s_dual8 = s_dual8value,
              r_primal9 = r_primal9value, s_dual9 = s_dual9value))
}
