library(Rcpp)
library(RcppArmadillo)
sourceCpp('~/Library/CloudStorage/GoogleDrive-96limtotoro@gmail.com/My Drive/research_multitask/multitask_supplement.cpp')

l2norm = function(x){
  sqrt(sum(x^2))
}

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

twoone_norm <- function(M){
  p = nrow(M)
  val = 0
  for(k in 1:p){
    val = val + sqrt(sum((M[k,]^2)))
  }
  return(val)
}
Frob_norm <- function(Mat){
  sqrt(sum(Mat^2))
}
Frob_norm_sq = function(Mat){
  sum(Mat^2)
}

rmse = function(x,y){
  sqrt(mean((x-y)^2))
}
sse = function(x,y){
  sum((x-y)^2)
}

sum_l2_norm_rows <- function(W) {
  sum(apply(W, 1, function(row) sqrt(sum(row^2))))
}

update_row_block <- function(Aj, Xj_norm2, lambda) {
  norm_Aj <- sqrt(sum(Aj^2))
  if (norm_Aj > lambda){
    return((1 - lambda / norm_Aj) * Aj / Xj_norm2)
  }else{
    return(rep(0, length(Aj)))
  }
}


MTL_admm = function(X, Y, lambda, 
                    rho = 1, max_iter = 1000, 
                    tol = 1e-5,
                    mu = 10, tau_inc = 2, tau_dec = 1.5,
                    verbose = FALSE) {
  
  start <- Sys.time()
  
  
  
  
  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(Y)
  
  B <- matrix(0, p, K)
  Gamma <- matrix(0, p, K)
  U <- matrix(0, p, K)  # scaled dual variable
  
  #if(verbose){print('making matrices')}
  XtX <- crossprod(X, X)
  XtY <- crossprod(X, Y)
  Ip <- diag(p)
  
  for (iter in 1:max_iter) {
    
    #if(verbose){print(iter)}

    Gamma_old <- Gamma
    
    # B-update (ridge-like step)
    B <- solve(1/(n*K) * XtX + rho * Ip, 1/(n*K) * XtY + rho * (Gamma - U))
    
    # Gamma-update (prox for \ell_2 norm: https://math.stackexchange.com/questions/2190885/proximal-operator-of-the-euclidean-norm-l-2-norm)
    Z <- B + U
    for (j in 1:p) {
      z_j <- Z[j, ]
      norm_z <- sqrt(sum(z_j^2))
      if (norm_z > lambda / rho) {
        Gamma[j, ] <- (1 - lambda / (rho * norm_z)) * z_j
      } else {
        Gamma[j, ] <- 0
      }
    }
    
    # Dual update
    U <- U + B - Gamma
    
    # Residuals
    primal_resid <- Frob_norm(B - Gamma)
    dual_resid <- rho * Frob_norm(Gamma - Gamma_old)
    
    if (verbose && iter %% 500 == 0) {
      cat(sprintf("Iter %d | Primal Resid: %.4f | Dual Resid: %.4f | Rho: %.3f\n",
                  iter, primal_resid, dual_resid, rho))
    }
    
    # Adaptive rho
    if (primal_resid > mu * dual_resid) {
      rho <- rho * tau_inc
      U <- U / tau_inc # change of variable required for the scaled variable
    } else if (dual_resid > mu * primal_resid) {
      rho <- rho / tau_dec
      U <- U * tau_dec # change of variable required for the scaled variable
    }
    
    
    # Check convergence
    if (iter > 2 & primal_resid < tol && dual_resid < tol) {
      if (verbose) cat(sprintf("Converged at iter %d\n", iter))
      break
    }
  }
  
  end  <- Sys.time()
  
  return(list(sol = Gamma,
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs"))))
}

admm_multitask <- function(X0, Y0, X1, Y1, 
                           X0tX0 = NULL, X0tY0 = NULL, X1tX1 = NULL, X1tY1 = NULL,
                           lambda0, lambda1, 
                           rho00 = 1, rho01 = 1, rho1 = 1, 
                           mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                           max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F){
  
  start <- Sys.time()
  
  
  
  
  
  
  N = nrow(X0) + nrow(X1)
  
  rho00value = rep(NA, max_iter)
  rho01value = rep(NA, max_iter)
  rho1value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  
  p <- ncol(X0)
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K)
  B1 <- matrix(0, p, K)
  Gamma00 <- matrix(0, p, K)
  Gamma01 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K)
  U0 <- matrix(0, p, K)
  U1 <- matrix(0, p, K)
  V1 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    X0tX0 = crossprod(X0,X0) ; X1tX1 = crossprod(X1,X1)
    X0tY0 = crossprod(X0,Y0) ; X1tY1 = crossprod(X1,Y1)
  }
  Ip = diag(p)
  
  for (iter in 1:max_iter) {
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma1_old = Gamma1
    
    #### Step 1: B0 update ####
    lhs = 1/(N*K) *X0tX0 + (rho00 + rho01)*Ip ; rhs = 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1)
    B0 <- solve(lhs, rhs)
    
    #### Step 2: B1 update ####
    lhs = 1/(N*K) *X1tX1 + rho1*Ip ; rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### Step 3: Gamma00 update ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)
    
    
    # #### Step 4: Gamma01 update ####
    # Gamma01 = update_Gamma01(B0,U1,Gamma1,lambda1,rho01,p,K)
    # #### Step 5: Gamma1 update ####
    # Gamma1 = update_Gamma1(B1,V1,Gamma01,lambda1,rho1,p,K)
    
    #### Step 4-5: exact joint update for (Gamma01, Gamma1) ####
    pair_update <- update_Gamma01_Gamma1_exact(B0, B1, U1, V1, lambda1, rho01, rho1, p, K)
    Gamma01 <- pair_update$Gamma01
    Gamma1  <- pair_update$Gamma1
    
    
    #### Step 6–8: Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    V1 <- V1 + B1 - Gamma1
    
    #### Convergence check: primal & dual ####
    r_primal00 = Frob_norm(B0 - Gamma00) ; s_dual00 = rho00 * Frob_norm(Gamma00 - Gamma00_old)
    r_primal01 = Frob_norm(B0 - Gamma01) ; s_dual01 = rho01 * Frob_norm(Gamma01 - Gamma01_old)
    r_primal1 = Frob_norm(B1 - Gamma1) ; s_dual1 = rho1 * Frob_norm(Gamma1 - Gamma1_old)
    
    
    
    if (verbose && iter %% 500 == 0) {
      cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f \n",
                  iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal1,base=10), log(s_dual1,base=10)))
    }
    
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
  
  end <- Sys.time()
  
  return(list(B0 = B0, B1 = B1, 
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma1 = Gamma1,
              U0 = U0, U1 = U1, V1 = V1, 
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
              rho00 = rho00value, rho01 = rho01value, rho1 = rho1value,
              r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal1 = r_primal1value, s_dual1 = s_dual1value))
}

admm_two_source <- function(X0, Y0, X1, Y1, X2, Y2,
                            X0tX0 = NULL, X0tY0 = NULL, X1tX1 = NULL, X1tY1 = NULL, X2tX2 = NULL, X2tY2 = NULL,
                            lambda0, lambda1, lambda2,
                            rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1, 
                            mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                            max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F, compute_objvalue = F){
  
  start <- Sys.time()
  
  N = nrow(X0) + nrow(X1) + nrow(X2)
  
  if(compute_objvalue == T){
    objvalue = rep(NA, max_iter)
  }
  
  
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
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K)
  B1 <- matrix(0, p, K)
  B2 <- matrix(0, p, K)
  Gamma00 <- matrix(0, p, K)
  Gamma01 <- matrix(0, p, K)
  Gamma02 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K)
  Gamma2 <- matrix(0, p, K)
  U0 <- matrix(0, p, K)
  U1 <- matrix(0, p, K)
  U2 <- matrix(0, p, K)
  V1 <- matrix(0, p, K)
  V2 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    X0tX0 <- crossprod(X0, X0)
    X1tX1 <- crossprod(X1, X1)
    X2tX2 <- crossprod(X2, X2)
    X0tY0 <- crossprod(X0, Y0)
    X1tY1 <- crossprod(X1, Y1)
    X2tY2 <- crossprod(X2, Y2)
  }
  Ip <- diag(p)
  
  
  for (iter in 1:max_iter){ # main iteration
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    
    #### Step 1: B0 update ####
    lhs <- 1/(N*K) *X0tX0 + (rho00 + rho01 + rho02) * Ip
    rhs <- 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2)
    B0 <- solve(lhs, rhs)
    
    #### Step 2: B1 update ####
    lhs = 1/(N*K) *X1tX1 + rho1*Ip
    rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### Step 3: B2 update ####
    lhs = 1/(N*K) *X2tX2 + rho2*Ip
    rhs = 1/(N*K) *X2tY2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
  
    #### Step 4: Gamma00 update ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)

    # #### Step 5: Gamma01 update ####
    # Gamma01 = update_Gamma01(B0,U1,Gamma1,lambda1,rho01,p,K)
    # #### Step 6: Gamma02 update ####
    # Gamma02 = update_Gamma01(B0,U2,Gamma2,lambda2,rho02,p,K)
    # #### Step 7: Gamma1 update ####
    # Gamma1 = update_Gamma1(B1,V1,Gamma01,lambda1,rho1,p,K)
    # #### Step 8: Gamma2 update ####
    # Gamma2 = update_Gamma1(B2,V2,Gamma02,lambda2,rho2,p,K)
    
    #### Step 5-6: exact joint update for (Gamma01, Gamma1) ####
    pair1 <- update_Gamma01_Gamma1_exact(B0, B1, U1, V1, lambda1, rho01, rho1, p, K)
    Gamma01 <- pair1$Gamma01
    Gamma1  <- pair1$Gamma1
    #### Step 7-8: exact joint update for (Gamma02, Gamma2) ####
    pair2 <- update_Gamma01_Gamma1_exact(B0, B2, U2, V2, lambda2, rho02, rho2, p, K)
    Gamma02 <- pair2$Gamma01
    Gamma2  <- pair2$Gamma1
    
    #### Step 9–13: Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    
    #### Convergence check: primal & dual ####
    r_primal00 = Frob_norm(B0 - Gamma00) ; s_dual00 = rho00 * Frob_norm(Gamma00 - Gamma00_old)
    r_primal01 = Frob_norm(B0 - Gamma01) ; s_dual01 = rho01 * Frob_norm(Gamma01 - Gamma01_old)
    r_primal02 = Frob_norm(B0 - Gamma02) ; s_dual02 = rho02 * Frob_norm(Gamma02 - Gamma02_old)
    
    r_primal1 = Frob_norm(B1 - Gamma1) ; s_dual1 = rho1 * Frob_norm(Gamma1 - Gamma1_old)
    r_primal2 = Frob_norm(B2 - Gamma2) ; s_dual2 = rho2 * Frob_norm(Gamma2 - Gamma2_old)
    
    if (verbose && iter %% 500 == 0) {
      cat(sprintf("Iter %d | Primal Resid00: %.2f | Dual Resid00: %.2f | Primal Resid01: %.2f | Dual Resid01: %.2f | Primal Resid02: %.2f | Dual Resid02: %.2f | Primal Resid1: %.2f | Dual Resid1: %.2f | Primal Resid2: %.2f | Dual Resid2: %.2f | \n",
                  iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10)))
    }
    
    if(compute_objvalue == T){
      objvalue[iter] = 1/(2*N) * (Frob_norm_sq(Y_target - X_target %*% B0) + Frob_norm_sq(Y_source1 - X_source1 %*% B1) + Frob_norm_sq(Y_source2 - X_source2 %*% B2)) +
        lambda0 * twoone_norm(B0) + lambda1 * twoone_norm(B0 - B1) + lambda2 * twoone_norm(B0 - B2)
    }
    
    
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
  
  end <- Sys.time()
  
  if(compute_objvalue == T){
    return(list(B0 = B0, B1 = B1, B2 = B2,
                Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma1 = Gamma1, Gamma2 = Gamma2,
                U0 = U0, U1 = U1, U2 = U2, V1 = V1, V2 = V2,
                iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
                rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho1 = rho1value, rho2 = rho2value,
                objvalue = objvalue,
                r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
                r_primal01 = r_primal01value, s_dual01 = s_dual01value,
                r_primal02 = r_primal02value, s_dual02 = s_dual02value,
                r_primal1 = r_primal1value, s_dual1 = s_dual1value,
                r_primal2 = r_primal2value, s_dual2 = s_dual2value))
  }else{
  return(list(B0 = B0, B1 = B1, B2 = B2,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma1 = Gamma1, Gamma2 = Gamma2,
              U0 = U0, U1 = U1, U2 = U2, V1 = V1, V2 = V2,
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho1 = rho1value, rho2 = rho2value,
              r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value))
  }
}

admm_three_source <- function(X0, Y0, X1, Y1, X2, Y2, X3, Y3,
                              X0tX0 = NULL, X0tY0 = NULL, X1tX1 = NULL, X1tY1 = NULL, X2tX2 = NULL, X2tY2 = NULL, X3tX3 = NULL, X3tY3 = NULL,
                              lambda0, lambda1, lambda2, lambda3,
                              rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho1 = 1, rho2 = 1, rho3 = 1,
                              mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                              max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F){
  
  start <- Sys.time()
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3)
  
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
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K)
  B1 <- matrix(0, p, K)
  B2 <- matrix(0, p, K)
  B3 <- matrix(0, p, K)
  Gamma00 <- matrix(0, p, K)
  Gamma01 <- matrix(0, p, K)
  Gamma02 <- matrix(0, p, K)
  Gamma03 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K)
  Gamma2 <- matrix(0, p, K)
  Gamma3 <- matrix(0, p, K)
  U0 <- matrix(0, p, K)
  U1 <- matrix(0, p, K)
  U2 <- matrix(0, p, K)
  U3 <- matrix(0, p, K)
  V1 <- matrix(0, p, K)
  V2 <- matrix(0, p, K)
  V3 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    X0tX0 <- crossprod(X0, X0)
    X1tX1 <- crossprod(X1, X1)
    X2tX2 <- crossprod(X2, X2)
    X3tX3 <- crossprod(X3, X3)
    X0tY0 <- crossprod(X0, Y0)
    X1tY1 <- crossprod(X1, Y1)
    X2tY2 <- crossprod(X2, Y2)
    X3tY3 <- crossprod(X3, Y3)
  }
  Ip <- diag(p)
  
  
  for (iter in 1:max_iter){ # main iteration
    
    Gamma00_old = Gamma00
    Gamma01_old = Gamma01
    Gamma02_old = Gamma02
    Gamma03_old = Gamma03
    Gamma1_old = Gamma1
    Gamma2_old = Gamma2
    Gamma3_old = Gamma3
    
    #### Step 1: B0 update ####
    lhs <- 1/(N*K) *X0tX0 + (rho00 + rho01 + rho02 + rho03) * Ip
    rhs <- 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3)
    B0 <- solve(lhs, rhs)
    
    #### Step 2: B1 update ####
    lhs = 1/(N*K) *X1tX1 + rho1*Ip
    rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### Step 3: B2 update ####
    lhs = 1/(N*K) *X2tX2 + rho2*Ip
    rhs = 1/(N*K) *X2tY2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### Step 4: B3 update ####
    lhs = 1/(N*K) *X3tX3 + rho3*Ip
    rhs = 1/(N*K) *X3tY3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### Step 5: Gamma00 update ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)
    # #### Step 6-8: Gamma01 - Gamma03 update ####
    # Gamma01 = update_Gamma01(B0,U1,Gamma1,lambda1,rho01,p,K)
    # Gamma02 = update_Gamma01(B0,U2,Gamma2,lambda2,rho02,p,K)
    # Gamma03 = update_Gamma01(B0,U3,Gamma3,lambda3,rho03,p,K)
    # 
    # #### Step 9-11: Gamma1 update ####
    # Gamma1 = update_Gamma1(B1,V1,Gamma01,lambda1,rho1,p,K)
    # Gamma2 = update_Gamma1(B2,V2,Gamma02,lambda2,rho2,p,K)
    # Gamma3 = update_Gamma1(B3,V3,Gamma03,lambda3,rho3,p,K)
    
    #### Step 6-7: exact joint update for (Gamma01, Gamma1) ####
    pair1 <- update_Gamma01_Gamma1_exact(B0, B1, U1, V1, lambda1, rho01, rho1, p, K)
    Gamma01 <- pair1$Gamma01
    Gamma1  <- pair1$Gamma1
    #### Step 8-9: exact joint update for (Gamma02, Gamma2) ####
    pair2 <- update_Gamma01_Gamma1_exact(B0, B2, U2, V2, lambda2, rho02, rho2, p, K)
    Gamma02 <- pair2$Gamma01
    Gamma2  <- pair2$Gamma1
    #### Step 10-11: exact joint update for (Gamma03, Gamma3) ####
    pair3 <- update_Gamma01_Gamma1_exact(B0, B3, U3, V3, lambda3, rho03, rho3, p, K)
    Gamma03 <- pair3$Gamma01
    Gamma3  <- pair3$Gamma1
    
    #### Step 12–18: Dual updates ####
    U0 <- U0 + B0 - Gamma00
    U1 <- U1 + B0 - Gamma01
    U2 <- U2 + B0 - Gamma02
    U3 <- U3 + B0 - Gamma03
    
    V1 <- V1 + B1 - Gamma1
    V2 <- V2 + B2 - Gamma2
    V3 <- V3 + B3 - Gamma3
    
    #### Convergence check: primal & dual ####
    r_primal00 = Frob_norm(B0 - Gamma00) ; s_dual00 = rho00 * Frob_norm(Gamma00 - Gamma00_old)
    r_primal01 = Frob_norm(B0 - Gamma01) ; s_dual01 = rho01 * Frob_norm(Gamma01 - Gamma01_old)
    r_primal02 = Frob_norm(B0 - Gamma02) ; s_dual02 = rho02 * Frob_norm(Gamma02 - Gamma02_old)
    r_primal03 = Frob_norm(B0 - Gamma03) ; s_dual03 = rho03 * Frob_norm(Gamma03 - Gamma03_old)
    
    r_primal1 = Frob_norm(B1 - Gamma1) ; s_dual1 = rho1 * Frob_norm(Gamma1 - Gamma1_old)
    r_primal2 = Frob_norm(B2 - Gamma2) ; s_dual2 = rho2 * Frob_norm(Gamma2 - Gamma2_old)
    r_primal3 = Frob_norm(B3 - Gamma3) ; s_dual3 = rho3 * Frob_norm(Gamma3 - Gamma3_old)
    
    if (verbose && iter %% 1000 == 0) {
      cat(sprintf("Iter %d | Primal Resid00: %.1f | Dual Resid00: %.1f | Primal Resid01: %.1f | Dual Resid01: %.1f | Primal Resid02: %.1f | Dual Resid02: %.1f | Primal Resid03: %.1f | Dual Resid03: %.1f | Primal Resid1: %.1f | Dual Resid1: %.1f | Primal Resid2: %.1f | Dual Resid2: %.1f | Primal Resid3: %.1f | Dual Resid3: %.1f | \n",
                  iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal03,base=10), log(s_dual03,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10), log(r_primal3,base=10), log(s_dual3,base=10)))
    }
    
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
  
  end <- Sys.time()
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, V1 = V1, V2 = V2, V3 = V3,
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, rho1 = rho1value, rho2 = rho2value, rho3 = rho3value,
              r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value))
}

admm_seven_source <- function(X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5, X6, Y6, X7, Y7,
                              X0tX0 = NULL, X0tY0 = NULL, X1tX1 = NULL, X1tY1 = NULL, X2tX2 = NULL, X2tY2 = NULL,
                              X3tX3 = NULL, X3tY3 = NULL, X4tX4 = NULL, X4tY4 = NULL, X5tX5 = NULL, X5tY5 = NULL, 
                              X6tX6 = NULL, X6tY6 = NULL, X7tX7 = NULL, X7tY7 = NULL, 
                              lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, 
                              rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1, rho07 = 1, 
                              rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1, rho7 = 1,
                              mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                              max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F){
  start <- Sys.time()
  
  N = nrow(X0) + nrow(X1) + nrow(X2) + nrow(X3) + nrow(X4) + nrow(X5) + nrow(X6) + nrow(X7)
  
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
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K)
  B1 <- matrix(0, p, K)
  B2 <- matrix(0, p, K)
  B3 <- matrix(0, p, K)
  B4 <- matrix(0, p, K)
  B5 <- matrix(0, p, K)
  B6 <- matrix(0, p, K)
  B7 <- matrix(0, p, K)
  
  Gamma00 <- matrix(0, p, K)
  Gamma01 <- matrix(0, p, K)
  Gamma02 <- matrix(0, p, K)
  Gamma03 <- matrix(0, p, K)
  Gamma04 <- matrix(0, p, K)
  Gamma05 <- matrix(0, p, K)
  Gamma06 <- matrix(0, p, K)
  Gamma07 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K)
  Gamma2 <- matrix(0, p, K)
  Gamma3 <- matrix(0, p, K)
  Gamma4 <- matrix(0, p, K)
  Gamma5 <- matrix(0, p, K)
  Gamma6 <- matrix(0, p, K)
  Gamma7 <- matrix(0, p, K)
  
  U0 <- matrix(0, p, K)
  U1 <- matrix(0, p, K)
  U2 <- matrix(0, p, K)
  U3 <- matrix(0, p, K)
  U4 <- matrix(0, p, K)
  U5 <- matrix(0, p, K)
  U6 <- matrix(0, p, K)
  U7 <- matrix(0, p, K)
  
  V1 <- matrix(0, p, K)
  V2 <- matrix(0, p, K)
  V3 <- matrix(0, p, K)
  V4 <- matrix(0, p, K)
  V5 <- matrix(0, p, K)
  V6 <- matrix(0, p, K)
  V7 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    X0tX0 <- crossprod(X0, X0)
    X1tX1 <- crossprod(X1, X1)
    X2tX2 <- crossprod(X2, X2)
    X3tX3 <- crossprod(X3, X3)
    X4tX4 <- crossprod(X4, X4) 
    X5tX5 <- crossprod(X5, X5)
    X6tX6 <- crossprod(X6, X6)
    X7tX7 <- crossprod(X7, X7)
    
    X0tY0 <- crossprod(X0, Y0)
    X1tY1 <- crossprod(X1, Y1)
    X2tY2 <- crossprod(X2, Y2)
    X3tY3 <- crossprod(X3, Y3)
    X4tY4 <- crossprod(X4, Y4)
    X5tY5 <- crossprod(X5, Y5)
    X6tY6 <- crossprod(X6, Y6)
    X7tY7 <- crossprod(X7, Y7)
  }
  Ip <- diag(p)
  
  
  for (iter in 1:max_iter){ # main iteration
    
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
    lhs <- 1/(N*K) *X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06 + rho07) * Ip
    rhs <- 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6) + rho07 * (Gamma07 - U7)
    B0 <- solve(lhs, rhs)
    
    #### Step 2: B1 update ####
    lhs = 1/(N*K) *X1tX1 + rho1*Ip
    rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### Step 3: B2 update ####
    lhs = 1/(N*K) *X2tX2 + rho2*Ip
    rhs = 1/(N*K) *X2tY2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### Step 4: B3 update ####
    lhs = 1/(N*K) *X3tX3 + rho3*Ip
    rhs = 1/(N*K) *X3tY3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### Step 5: B4 update ####
    lhs = 1/(N*K) *X4tX4 + rho4*Ip
    rhs = 1/(N*K) *X4tY4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### Step 6: B5 update ####
    lhs = 1/(N*K) *X5tX5 + rho5*Ip
    rhs = 1/(N*K) *X5tY5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    #### B6 update ####
    lhs = 1/(N*K) *X6tX6 + rho6*Ip
    rhs = 1/(N*K) *X6tY6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    #### B7 update ####
    lhs = 1/(N*K) *X7tX7 + rho7*Ip
    rhs = 1/(N*K) *X7tY7 + rho7*(Gamma7 - V7)
    B7 <- solve(lhs, rhs)
    
    ####  ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)
    # ####  ####
    # Gamma01 = update_Gamma01(B0,U1,Gamma1,lambda1,rho01,p,K)
    # Gamma02 = update_Gamma01(B0,U2,Gamma2,lambda2,rho02,p,K)
    # Gamma03 = update_Gamma01(B0,U3,Gamma3,lambda3,rho03,p,K)
    # Gamma04 = update_Gamma01(B0,U4,Gamma4,lambda4,rho04,p,K)
    # Gamma05 = update_Gamma01(B0,U5,Gamma5,lambda5,rho05,p,K)
    # Gamma06 = update_Gamma01(B0,U6,Gamma6,lambda6,rho06,p,K)
    # Gamma07 = update_Gamma01(B0,U7,Gamma7,lambda7,rho07,p,K)
    # 
    # ####  ####
    # Gamma1 = update_Gamma1(B1,V1,Gamma01,lambda1,rho1,p,K)
    # Gamma2 = update_Gamma1(B2,V2,Gamma02,lambda2,rho2,p,K)
    # Gamma3 = update_Gamma1(B3,V3,Gamma03,lambda3,rho3,p,K)
    # Gamma4 = update_Gamma1(B4,V4,Gamma04,lambda4,rho4,p,K)
    # Gamma5 = update_Gamma1(B5,V5,Gamma05,lambda5,rho5,p,K)
    # Gamma6 = update_Gamma1(B6,V6,Gamma06,lambda6,rho6,p,K)
    # Gamma7 = update_Gamma1(B7,V7,Gamma07,lambda7,rho7,p,K)
    ####  ####
    pair1 <- update_Gamma01_Gamma1_exact(B0, B1, U1, V1, lambda1, rho01, rho1, p, K)
    Gamma01 <- pair1$Gamma01
    Gamma1  <- pair1$Gamma1
    
    pair2 <- update_Gamma01_Gamma1_exact(B0, B2, U2, V2, lambda2, rho02, rho2, p, K)
    Gamma02 <- pair2$Gamma01
    Gamma2  <- pair2$Gamma1
    
    pair3 <- update_Gamma01_Gamma1_exact(B0, B3, U3, V3, lambda3, rho03, rho3, p, K)
    Gamma03 <- pair3$Gamma01
    Gamma3  <- pair3$Gamma1
    
    pair4 <- update_Gamma01_Gamma1_exact(B0, B4, U4, V4, lambda4, rho04, rho4, p, K)
    Gamma04 <- pair4$Gamma01
    Gamma4  <- pair4$Gamma1
    
    pair5 <- update_Gamma01_Gamma1_exact(B0, B5, U5, V5, lambda5, rho05, rho5, p, K)
    Gamma05 <- pair5$Gamma01
    Gamma5  <- pair5$Gamma1
    
    pair6 <- update_Gamma01_Gamma1_exact(B0, B6, U6, V6, lambda6, rho06, rho6, p, K)
    Gamma06 <- pair6$Gamma01
    Gamma6  <- pair6$Gamma1
    
    pair7 <- update_Gamma01_Gamma1_exact(B0, B7, U7, V7, lambda7, rho07, rho7, p, K)
    Gamma07 <- pair7$Gamma01
    Gamma7  <- pair7$Gamma1
    
    ####  Dual updates ####
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
    r_primal00 = Frob_norm(B0 - Gamma00) ; s_dual00 = rho00 * Frob_norm(Gamma00 - Gamma00_old)
    r_primal01 = Frob_norm(B0 - Gamma01) ; s_dual01 = rho01 * Frob_norm(Gamma01 - Gamma01_old)
    r_primal02 = Frob_norm(B0 - Gamma02) ; s_dual02 = rho02 * Frob_norm(Gamma02 - Gamma02_old)
    r_primal03 = Frob_norm(B0 - Gamma03) ; s_dual03 = rho03 * Frob_norm(Gamma03 - Gamma03_old)
    r_primal04 = Frob_norm(B0 - Gamma04) ; s_dual04 = rho04 * Frob_norm(Gamma04 - Gamma04_old)
    r_primal05 = Frob_norm(B0 - Gamma05) ; s_dual05 = rho05 * Frob_norm(Gamma05 - Gamma05_old)
    r_primal06 = Frob_norm(B0 - Gamma06) ; s_dual06 = rho06 * Frob_norm(Gamma06 - Gamma06_old)
    r_primal07 = Frob_norm(B0 - Gamma07) ; s_dual07 = rho07 * Frob_norm(Gamma07 - Gamma07_old)
    
    r_primal1 = Frob_norm(B1 - Gamma1) ; s_dual1 = rho1 * Frob_norm(Gamma1 - Gamma1_old)
    r_primal2 = Frob_norm(B2 - Gamma2) ; s_dual2 = rho2 * Frob_norm(Gamma2 - Gamma2_old)
    r_primal3 = Frob_norm(B3 - Gamma3) ; s_dual3 = rho3 * Frob_norm(Gamma3 - Gamma3_old)
    r_primal4 = Frob_norm(B4 - Gamma4) ; s_dual4 = rho4 * Frob_norm(Gamma4 - Gamma4_old)
    r_primal5 = Frob_norm(B5 - Gamma5) ; s_dual5 = rho5 * Frob_norm(Gamma5 - Gamma5_old)
    r_primal6 = Frob_norm(B6 - Gamma6) ; s_dual6 = rho6 * Frob_norm(Gamma6 - Gamma6_old)
    r_primal7 = Frob_norm(B7 - Gamma7) ; s_dual7 = rho7 * Frob_norm(Gamma7 - Gamma7_old)
    
    
    # if (verbose && iter %% 1000 == 0) {
    #   cat(sprintf("Iter %d | Primal Resid00: %.1f | Dual Resid00: %.1f | Primal Resid01: %.1f | Dual Resid01: %.1f | Primal Resid02: %.1f | Dual Resid02: %.1f | Primal Resid03: %.1f | Dual Resid03: %.1f | Primal Resid1: %.1f | Dual Resid1: %.1f | Primal Resid2: %.1f | Dual Resid2: %.1f | Primal Resid3: %.1f | Dual Resid3: %.1f | \n",
    #               iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal03,base=10), log(s_dual03,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10), log(r_primal3,base=10), log(s_dual3,base=10)))
    # }
    
    
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
  
  end <- Sys.time()
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6, B7 = B7,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06, Gamma07 = Gamma07,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6, Gamma7 = Gamma7,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6, U7 = U7,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6, V7 = V7,
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
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

# Additional functions from functions_joint_cluster.R

admm_eleven_source <- function(X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5, X6, Y6, X7, Y7, X8, Y8, X9, Y9, X10, Y10, X11, Y11,
                               X0tX0 = NULL, X0tY0 = NULL,
                               X1tX1 = NULL, X1tY1 = NULL,
                               X2tX2 = NULL, X2tY2 = NULL,
                               X3tX3 = NULL, X3tY3 = NULL,
                               X4tX4 = NULL, X4tY4 = NULL,
                               X5tX5 = NULL, X5tY5 = NULL,
                               X6tX6 = NULL, X6tY6 = NULL,
                               X7tX7 = NULL, X7tY7 = NULL,
                               X8tX8 = NULL, X8tY8 = NULL,
                               X9tX9 = NULL, X9tY9 = NULL,
                               X10tX10 = NULL, X10tY10 = NULL,
                               X11tX11 = NULL, X11tY11 = NULL,
                               lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8, lambda9, lambda10, lambda11,
                               rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1, rho07 = 1, rho08 = 1, rho09 = 1, rho010 = 1, rho011 = 1,
                               rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1, rho7 = 1, rho8 = 1, rho9 = 1, rho10 = 1, rho11 = 1,
                               mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                               max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F){
  start <- Sys.time()
  
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
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K)
  B1 <- matrix(0, p, K)
  B2 <- matrix(0, p, K)
  B3 <- matrix(0, p, K)
  B4 <- matrix(0, p, K)
  B5 <- matrix(0, p, K)
  B6 <- matrix(0, p, K)
  B7 <- matrix(0, p, K)
  B8 <- matrix(0, p, K)
  B9 <- matrix(0, p, K)
  B10 <- matrix(0, p, K)
  B11 <- matrix(0, p, K)
  
  Gamma00 <- matrix(0, p, K)
  Gamma01 <- matrix(0, p, K)
  Gamma02 <- matrix(0, p, K)
  Gamma03 <- matrix(0, p, K)
  Gamma04 <- matrix(0, p, K)
  Gamma05 <- matrix(0, p, K)
  Gamma06 <- matrix(0, p, K)
  Gamma07 <- matrix(0, p, K)
  Gamma08 <- matrix(0, p, K)
  Gamma09 <- matrix(0, p, K)
  Gamma010 <- matrix(0, p, K)
  Gamma011 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K)
  Gamma2 <- matrix(0, p, K)
  Gamma3 <- matrix(0, p, K)
  Gamma4 <- matrix(0, p, K)
  Gamma5 <- matrix(0, p, K)
  Gamma6 <- matrix(0, p, K)
  Gamma7 <- matrix(0, p, K)
  Gamma8 <- matrix(0, p, K)
  Gamma9 <- matrix(0, p, K)
  Gamma10 <- matrix(0, p, K)
  Gamma11 <- matrix(0, p, K)
  
  U0 <- matrix(0, p, K)
  U1 <- matrix(0, p, K)
  U2 <- matrix(0, p, K)
  U3 <- matrix(0, p, K)
  U4 <- matrix(0, p, K)
  U5 <- matrix(0, p, K)
  U6 <- matrix(0, p, K)
  U7 <- matrix(0, p, K)
  U8 <- matrix(0, p, K)
  U9 <- matrix(0, p, K)
  U10 <- matrix(0, p, K)
  U11 <- matrix(0, p, K)
  
  V1 <- matrix(0, p, K)
  V2 <- matrix(0, p, K)
  V3 <- matrix(0, p, K)
  V4 <- matrix(0, p, K)
  V5 <- matrix(0, p, K)
  V6 <- matrix(0, p, K)
  V7 <- matrix(0, p, K)
  V8 <- matrix(0, p, K)
  V9 <- matrix(0, p, K)
  V10 <- matrix(0, p, K)
  V11 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    X0tX0 <- crossprod(X0, X0)
    X1tX1 <- crossprod(X1, X1)
    X2tX2 <- crossprod(X2, X2)
    X3tX3 <- crossprod(X3, X3)
    X4tX4 <- crossprod(X4, X4)
    X5tX5 <- crossprod(X5, X5)
    X6tX6 <- crossprod(X6, X6)
    X7tX7 <- crossprod(X7, X7)
    X8tX8 <- crossprod(X8, X8)
    X9tX9 <- crossprod(X9, X9)
    X10tX10 <- crossprod(X10, X10)
    X11tX11 <- crossprod(X11, X11)
    
    X0tY0 <- crossprod(X0, Y0)
    X1tY1 <- crossprod(X1, Y1)
    X2tY2 <- crossprod(X2, Y2)
    X3tY3 <- crossprod(X3, Y3)
    X4tY4 <- crossprod(X4, Y4)
    X5tY5 <- crossprod(X5, Y5)
    X6tY6 <- crossprod(X6, Y6)
    X7tY7 <- crossprod(X7, Y7)
    X8tY8 <- crossprod(X8, Y8)
    X9tY9 <- crossprod(X9, Y9)
    X10tY10 <- crossprod(X10, Y10)
    X11tY11 <- crossprod(X11, Y11)
  }
  Ip <- diag(p)
  
  for (iter in 1:max_iter){ # main iteration
    
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
    lhs <- 1/(N*K) *X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06 + rho07 + rho08 + rho09 + rho010 + rho011) * Ip
    rhs <- 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6) + rho07 * (Gamma07 - U7) + rho08 * (Gamma08 - U8) + rho09 * (Gamma09 - U9) + rho010 * (Gamma010 - U10) + rho011 * (Gamma011 - U11)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/(N*K) *X1tX1 + rho1*Ip
    rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/(N*K) *X2tX2 + rho2*Ip
    rhs = 1/(N*K) *X2tY2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/(N*K) *X3tX3 + rho3*Ip
    rhs = 1/(N*K) *X3tY3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/(N*K) *X4tX4 + rho4*Ip
    rhs = 1/(N*K) *X4tY4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### B5 update ####
    lhs = 1/(N*K) *X5tX5 + rho5*Ip
    rhs = 1/(N*K) *X5tY5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    #### B6 update ####
    lhs = 1/(N*K) *X6tX6 + rho6*Ip
    rhs = 1/(N*K) *X6tY6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    
    #### B7 update ####
    lhs = 1/(N*K) *X7tX7 + rho7*Ip
    rhs = 1/(N*K) *X7tY7 + rho7*(Gamma7 - V7)
    B7 <- solve(lhs, rhs)
    
    #### B8 update ####
    lhs = 1/(N*K) *X8tX8 + rho8*Ip
    rhs = 1/(N*K) *X8tY8 + rho8*(Gamma8 - V8)
    B8 <- solve(lhs, rhs)
    
    #### B9 update ####
    lhs = 1/(N*K) *X9tX9 + rho9*Ip
    rhs = 1/(N*K) *X9tY9 + rho9*(Gamma9 - V9)
    B9 <- solve(lhs, rhs)
    
    #### B10 update ####
    lhs = 1/(N*K) *X10tX10 + rho10*Ip
    rhs = 1/(N*K) *X10tY10 + rho10*(Gamma10 - V10)
    B10 <- solve(lhs, rhs)
    
    #### B11 update ####
    lhs = 1/(N*K) *X11tX11 + rho11*Ip
    rhs = 1/(N*K) *X11tY11 + rho11*(Gamma11 - V11)
    B11 <- solve(lhs, rhs)
    
    ####  ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)
    
    pair1 <- update_Gamma01_Gamma1_exact(B0, B1, U1, V1, lambda1, rho01, rho1, p, K)
    Gamma01 <- pair1$Gamma01
    Gamma1  <- pair1$Gamma1
    
    pair2 <- update_Gamma01_Gamma1_exact(B0, B2, U2, V2, lambda2, rho02, rho2, p, K)
    Gamma02 <- pair2$Gamma01
    Gamma2  <- pair2$Gamma1
    
    pair3 <- update_Gamma01_Gamma1_exact(B0, B3, U3, V3, lambda3, rho03, rho3, p, K)
    Gamma03 <- pair3$Gamma01
    Gamma3  <- pair3$Gamma1
    
    pair4 <- update_Gamma01_Gamma1_exact(B0, B4, U4, V4, lambda4, rho04, rho4, p, K)
    Gamma04 <- pair4$Gamma01
    Gamma4  <- pair4$Gamma1
    
    pair5 <- update_Gamma01_Gamma1_exact(B0, B5, U5, V5, lambda5, rho05, rho5, p, K)
    Gamma05 <- pair5$Gamma01
    Gamma5  <- pair5$Gamma1
    
    pair6 <- update_Gamma01_Gamma1_exact(B0, B6, U6, V6, lambda6, rho06, rho6, p, K)
    Gamma06 <- pair6$Gamma01
    Gamma6  <- pair6$Gamma1
    
    pair7 <- update_Gamma01_Gamma1_exact(B0, B7, U7, V7, lambda7, rho07, rho7, p, K)
    Gamma07 <- pair7$Gamma01
    Gamma7  <- pair7$Gamma1
    
    pair8 <- update_Gamma01_Gamma1_exact(B0, B8, U8, V8, lambda8, rho08, rho8, p, K)
    Gamma08 <- pair8$Gamma01
    Gamma8  <- pair8$Gamma1
    
    pair9 <- update_Gamma01_Gamma1_exact(B0, B9, U9, V9, lambda9, rho09, rho9, p, K)
    Gamma09 <- pair9$Gamma01
    Gamma9  <- pair9$Gamma1
    
    pair10 <- update_Gamma01_Gamma1_exact(B0, B10, U10, V10, lambda10, rho010, rho10, p, K)
    Gamma010 <- pair10$Gamma01
    Gamma10  <- pair10$Gamma1
    
    pair11 <- update_Gamma01_Gamma1_exact(B0, B11, U11, V11, lambda11, rho011, rho11, p, K)
    Gamma011 <- pair11$Gamma01
    Gamma11  <- pair11$Gamma1
    
    ####  Dual updates ####
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
    r_primal00 = Frob_norm(B0 - Gamma00) ; s_dual00 = rho00 * Frob_norm(Gamma00 - Gamma00_old)
    r_primal01 = Frob_norm(B0 - Gamma01) ; s_dual01 = rho01 * Frob_norm(Gamma01 - Gamma01_old)
    r_primal02 = Frob_norm(B0 - Gamma02) ; s_dual02 = rho02 * Frob_norm(Gamma02 - Gamma02_old)
    r_primal03 = Frob_norm(B0 - Gamma03) ; s_dual03 = rho03 * Frob_norm(Gamma03 - Gamma03_old)
    r_primal04 = Frob_norm(B0 - Gamma04) ; s_dual04 = rho04 * Frob_norm(Gamma04 - Gamma04_old)
    r_primal05 = Frob_norm(B0 - Gamma05) ; s_dual05 = rho05 * Frob_norm(Gamma05 - Gamma05_old)
    r_primal06 = Frob_norm(B0 - Gamma06) ; s_dual06 = rho06 * Frob_norm(Gamma06 - Gamma06_old)
    r_primal07 = Frob_norm(B0 - Gamma07) ; s_dual07 = rho07 * Frob_norm(Gamma07 - Gamma07_old)
    r_primal08 = Frob_norm(B0 - Gamma08) ; s_dual08 = rho08 * Frob_norm(Gamma08 - Gamma08_old)
    r_primal09 = Frob_norm(B0 - Gamma09) ; s_dual09 = rho09 * Frob_norm(Gamma09 - Gamma09_old)
    r_primal010 = Frob_norm(B0 - Gamma010) ; s_dual010 = rho010 * Frob_norm(Gamma010 - Gamma010_old)
    r_primal011 = Frob_norm(B0 - Gamma011) ; s_dual011 = rho011 * Frob_norm(Gamma011 - Gamma011_old)
    
    r_primal1 = Frob_norm(B1 - Gamma1) ; s_dual1 = rho1 * Frob_norm(Gamma1 - Gamma1_old)
    r_primal2 = Frob_norm(B2 - Gamma2) ; s_dual2 = rho2 * Frob_norm(Gamma2 - Gamma2_old)
    r_primal3 = Frob_norm(B3 - Gamma3) ; s_dual3 = rho3 * Frob_norm(Gamma3 - Gamma3_old)
    r_primal4 = Frob_norm(B4 - Gamma4) ; s_dual4 = rho4 * Frob_norm(Gamma4 - Gamma4_old)
    r_primal5 = Frob_norm(B5 - Gamma5) ; s_dual5 = rho5 * Frob_norm(Gamma5 - Gamma5_old)
    r_primal6 = Frob_norm(B6 - Gamma6) ; s_dual6 = rho6 * Frob_norm(Gamma6 - Gamma6_old)
    r_primal7 = Frob_norm(B7 - Gamma7) ; s_dual7 = rho7 * Frob_norm(Gamma7 - Gamma7_old)
    r_primal8 = Frob_norm(B8 - Gamma8) ; s_dual8 = rho8 * Frob_norm(Gamma8 - Gamma8_old)
    r_primal9 = Frob_norm(B9 - Gamma9) ; s_dual9 = rho9 * Frob_norm(Gamma9 - Gamma9_old)
    r_primal10 = Frob_norm(B10 - Gamma10) ; s_dual10 = rho10 * Frob_norm(Gamma10 - Gamma10_old)
    r_primal11 = Frob_norm(B11 - Gamma11) ; s_dual11 = rho11 * Frob_norm(Gamma11 - Gamma11_old)
    
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
  
  end <- Sys.time()
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6, B7 = B7, B8 = B8, B9 = B9, B10 = B10, B11 = B11,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06, Gamma07 = Gamma07, Gamma08 = Gamma08, Gamma09 = Gamma09, Gamma010 = Gamma010, Gamma011 = Gamma011,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6, Gamma7 = Gamma7, Gamma8 = Gamma8, Gamma9 = Gamma9, Gamma10 = Gamma10, Gamma11 = Gamma11,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6, U7 = U7, U8 = U8, U9 = U9, U10 = U10, U11 = U11,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6, V7 = V7, V8 = V8, V9 = V9, V10 = V10, V11 = V11,
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
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

admm_four_source <- function(X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4,
                             X0tX0 = NULL, X0tY0 = NULL,
                             X1tX1 = NULL, X1tY1 = NULL,
                             X2tX2 = NULL, X2tY2 = NULL,
                             X3tX3 = NULL, X3tY3 = NULL,
                             X4tX4 = NULL, X4tY4 = NULL,
                             lambda0, lambda1, lambda2, lambda3, lambda4,
                             rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1,
                             rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1,
                             mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                             max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F){
  start <- Sys.time()
  
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
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K)
  B1 <- matrix(0, p, K)
  B2 <- matrix(0, p, K)
  B3 <- matrix(0, p, K)
  B4 <- matrix(0, p, K)
  
  Gamma00 <- matrix(0, p, K)
  Gamma01 <- matrix(0, p, K)
  Gamma02 <- matrix(0, p, K)
  Gamma03 <- matrix(0, p, K)
  Gamma04 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K)
  Gamma2 <- matrix(0, p, K)
  Gamma3 <- matrix(0, p, K)
  Gamma4 <- matrix(0, p, K)
  
  U0 <- matrix(0, p, K)
  U1 <- matrix(0, p, K)
  U2 <- matrix(0, p, K)
  U3 <- matrix(0, p, K)
  U4 <- matrix(0, p, K)
  
  V1 <- matrix(0, p, K)
  V2 <- matrix(0, p, K)
  V3 <- matrix(0, p, K)
  V4 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    X0tX0 <- crossprod(X0, X0)
    X1tX1 <- crossprod(X1, X1)
    X2tX2 <- crossprod(X2, X2)
    X3tX3 <- crossprod(X3, X3)
    X4tX4 <- crossprod(X4, X4)
    
    X0tY0 <- crossprod(X0, Y0)
    X1tY1 <- crossprod(X1, Y1)
    X2tY2 <- crossprod(X2, Y2)
    X3tY3 <- crossprod(X3, Y3)
    X4tY4 <- crossprod(X4, Y4)
  }
  Ip <- diag(p)
  
  for (iter in 1:max_iter){ # main iteration
    
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
    lhs <- 1/(N*K) *X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04) * Ip
    rhs <- 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/(N*K) *X1tX1 + rho1*Ip
    rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/(N*K) *X2tX2 + rho2*Ip
    rhs = 1/(N*K) *X2tY2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/(N*K) *X3tX3 + rho3*Ip
    rhs = 1/(N*K) *X3tY3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/(N*K) *X4tX4 + rho4*Ip
    rhs = 1/(N*K) *X4tY4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    ####  ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)
    
    pair1 <- update_Gamma01_Gamma1_exact(B0, B1, U1, V1, lambda1, rho01, rho1, p, K)
    Gamma01 <- pair1$Gamma01
    Gamma1  <- pair1$Gamma1
    
    pair2 <- update_Gamma01_Gamma1_exact(B0, B2, U2, V2, lambda2, rho02, rho2, p, K)
    Gamma02 <- pair2$Gamma01
    Gamma2  <- pair2$Gamma1
    
    pair3 <- update_Gamma01_Gamma1_exact(B0, B3, U3, V3, lambda3, rho03, rho3, p, K)
    Gamma03 <- pair3$Gamma01
    Gamma3  <- pair3$Gamma1
    
    pair4 <- update_Gamma01_Gamma1_exact(B0, B4, U4, V4, lambda4, rho04, rho4, p, K)
    Gamma04 <- pair4$Gamma01
    Gamma4  <- pair4$Gamma1
    
    ####  Dual updates ####
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
    r_primal00 = Frob_norm(B0 - Gamma00) ; s_dual00 = rho00 * Frob_norm(Gamma00 - Gamma00_old)
    r_primal01 = Frob_norm(B0 - Gamma01) ; s_dual01 = rho01 * Frob_norm(Gamma01 - Gamma01_old)
    r_primal02 = Frob_norm(B0 - Gamma02) ; s_dual02 = rho02 * Frob_norm(Gamma02 - Gamma02_old)
    r_primal03 = Frob_norm(B0 - Gamma03) ; s_dual03 = rho03 * Frob_norm(Gamma03 - Gamma03_old)
    r_primal04 = Frob_norm(B0 - Gamma04) ; s_dual04 = rho04 * Frob_norm(Gamma04 - Gamma04_old)
    
    r_primal1 = Frob_norm(B1 - Gamma1) ; s_dual1 = rho1 * Frob_norm(Gamma1 - Gamma1_old)
    r_primal2 = Frob_norm(B2 - Gamma2) ; s_dual2 = rho2 * Frob_norm(Gamma2 - Gamma2_old)
    r_primal3 = Frob_norm(B3 - Gamma3) ; s_dual3 = rho3 * Frob_norm(Gamma3 - Gamma3_old)
    r_primal4 = Frob_norm(B4 - Gamma4) ; s_dual4 = rho4 * Frob_norm(Gamma4 - Gamma4_old)
    
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
  
  end <- Sys.time()
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4,
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
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

admm_thirteen_source <- function(X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5, X6, Y6, X7, Y7, X8, Y8, X9, Y9, X10, Y10, X11, Y11, X12, Y12, X13, Y13,
                                 X0tX0 = NULL, X0tY0 = NULL,
                                 X1tX1 = NULL, X1tY1 = NULL,
                                 X2tX2 = NULL, X2tY2 = NULL,
                                 X3tX3 = NULL, X3tY3 = NULL,
                                 X4tX4 = NULL, X4tY4 = NULL,
                                 X5tX5 = NULL, X5tY5 = NULL,
                                 X6tX6 = NULL, X6tY6 = NULL,
                                 X7tX7 = NULL, X7tY7 = NULL,
                                 X8tX8 = NULL, X8tY8 = NULL,
                                 X9tX9 = NULL, X9tY9 = NULL,
                                 X10tX10 = NULL, X10tY10 = NULL,
                                 X11tX11 = NULL, X11tY11 = NULL,
                                 X12tX12 = NULL, X12tY12 = NULL,
                                 X13tX13 = NULL, X13tY13 = NULL,
                                 lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8, lambda9, lambda10, lambda11, lambda12, lambda13,
                                 rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1, rho07 = 1, rho08 = 1, rho09 = 1, rho010 = 1, rho011 = 1, rho012 = 1, rho013 = 1,
                                 rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1, rho7 = 1, rho8 = 1, rho9 = 1, rho10 = 1, rho11 = 1, rho12 = 1, rho13 = 1,
                                 mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                                 max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F){
  start <- Sys.time()
  
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
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K)
  B1 <- matrix(0, p, K)
  B2 <- matrix(0, p, K)
  B3 <- matrix(0, p, K)
  B4 <- matrix(0, p, K)
  B5 <- matrix(0, p, K)
  B6 <- matrix(0, p, K)
  B7 <- matrix(0, p, K)
  B8 <- matrix(0, p, K)
  B9 <- matrix(0, p, K)
  B10 <- matrix(0, p, K)
  B11 <- matrix(0, p, K)
  B12 <- matrix(0, p, K)
  B13 <- matrix(0, p, K)
  
  Gamma00 <- matrix(0, p, K)
  Gamma01 <- matrix(0, p, K)
  Gamma02 <- matrix(0, p, K)
  Gamma03 <- matrix(0, p, K)
  Gamma04 <- matrix(0, p, K)
  Gamma05 <- matrix(0, p, K)
  Gamma06 <- matrix(0, p, K)
  Gamma07 <- matrix(0, p, K)
  Gamma08 <- matrix(0, p, K)
  Gamma09 <- matrix(0, p, K)
  Gamma010 <- matrix(0, p, K)
  Gamma011 <- matrix(0, p, K)
  Gamma012 <- matrix(0, p, K)
  Gamma013 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K)
  Gamma2 <- matrix(0, p, K)
  Gamma3 <- matrix(0, p, K)
  Gamma4 <- matrix(0, p, K)
  Gamma5 <- matrix(0, p, K)
  Gamma6 <- matrix(0, p, K)
  Gamma7 <- matrix(0, p, K)
  Gamma8 <- matrix(0, p, K)
  Gamma9 <- matrix(0, p, K)
  Gamma10 <- matrix(0, p, K)
  Gamma11 <- matrix(0, p, K)
  Gamma12 <- matrix(0, p, K)
  Gamma13 <- matrix(0, p, K)
  
  U0 <- matrix(0, p, K)
  U1 <- matrix(0, p, K)
  U2 <- matrix(0, p, K)
  U3 <- matrix(0, p, K)
  U4 <- matrix(0, p, K)
  U5 <- matrix(0, p, K)
  U6 <- matrix(0, p, K)
  U7 <- matrix(0, p, K)
  U8 <- matrix(0, p, K)
  U9 <- matrix(0, p, K)
  U10 <- matrix(0, p, K)
  U11 <- matrix(0, p, K)
  U12 <- matrix(0, p, K)
  U13 <- matrix(0, p, K)
  
  V1 <- matrix(0, p, K)
  V2 <- matrix(0, p, K)
  V3 <- matrix(0, p, K)
  V4 <- matrix(0, p, K)
  V5 <- matrix(0, p, K)
  V6 <- matrix(0, p, K)
  V7 <- matrix(0, p, K)
  V8 <- matrix(0, p, K)
  V9 <- matrix(0, p, K)
  V10 <- matrix(0, p, K)
  V11 <- matrix(0, p, K)
  V12 <- matrix(0, p, K)
  V13 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    X0tX0 <- crossprod(X0, X0)
    X1tX1 <- crossprod(X1, X1)
    X2tX2 <- crossprod(X2, X2)
    X3tX3 <- crossprod(X3, X3)
    X4tX4 <- crossprod(X4, X4)
    X5tX5 <- crossprod(X5, X5)
    X6tX6 <- crossprod(X6, X6)
    X7tX7 <- crossprod(X7, X7)
    X8tX8 <- crossprod(X8, X8)
    X9tX9 <- crossprod(X9, X9)
    X10tX10 <- crossprod(X10, X10)
    X11tX11 <- crossprod(X11, X11)
    X12tX12 <- crossprod(X12, X12)
    X13tX13 <- crossprod(X13, X13)
    
    X0tY0 <- crossprod(X0, Y0)
    X1tY1 <- crossprod(X1, Y1)
    X2tY2 <- crossprod(X2, Y2)
    X3tY3 <- crossprod(X3, Y3)
    X4tY4 <- crossprod(X4, Y4)
    X5tY5 <- crossprod(X5, Y5)
    X6tY6 <- crossprod(X6, Y6)
    X7tY7 <- crossprod(X7, Y7)
    X8tY8 <- crossprod(X8, Y8)
    X9tY9 <- crossprod(X9, Y9)
    X10tY10 <- crossprod(X10, Y10)
    X11tY11 <- crossprod(X11, Y11)
    X12tY12 <- crossprod(X12, Y12)
    X13tY13 <- crossprod(X13, Y13)
  }
  Ip <- diag(p)
  
  for (iter in 1:max_iter){ # main iteration
    
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
    lhs <- 1/(N*K) *X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06 + rho07 + rho08 + rho09 + rho010 + rho011 + rho012 + rho013) * Ip
    rhs <- 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6) + rho07 * (Gamma07 - U7) + rho08 * (Gamma08 - U8) + rho09 * (Gamma09 - U9) + rho010 * (Gamma010 - U10) + rho011 * (Gamma011 - U11) + rho012 * (Gamma012 - U12) + rho013 * (Gamma013 - U13)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/(N*K) *X1tX1 + rho1*Ip
    rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/(N*K) *X2tX2 + rho2*Ip
    rhs = 1/(N*K) *X2tY2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/(N*K) *X3tX3 + rho3*Ip
    rhs = 1/(N*K) *X3tY3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/(N*K) *X4tX4 + rho4*Ip
    rhs = 1/(N*K) *X4tY4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### B5 update ####
    lhs = 1/(N*K) *X5tX5 + rho5*Ip
    rhs = 1/(N*K) *X5tY5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    #### B6 update ####
    lhs = 1/(N*K) *X6tX6 + rho6*Ip
    rhs = 1/(N*K) *X6tY6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    
    #### B7 update ####
    lhs = 1/(N*K) *X7tX7 + rho7*Ip
    rhs = 1/(N*K) *X7tY7 + rho7*(Gamma7 - V7)
    B7 <- solve(lhs, rhs)
    
    #### B8 update ####
    lhs = 1/(N*K) *X8tX8 + rho8*Ip
    rhs = 1/(N*K) *X8tY8 + rho8*(Gamma8 - V8)
    B8 <- solve(lhs, rhs)
    
    #### B9 update ####
    lhs = 1/(N*K) *X9tX9 + rho9*Ip
    rhs = 1/(N*K) *X9tY9 + rho9*(Gamma9 - V9)
    B9 <- solve(lhs, rhs)
    
    #### B10 update ####
    lhs = 1/(N*K) *X10tX10 + rho10*Ip
    rhs = 1/(N*K) *X10tY10 + rho10*(Gamma10 - V10)
    B10 <- solve(lhs, rhs)
    
    #### B11 update ####
    lhs = 1/(N*K) *X11tX11 + rho11*Ip
    rhs = 1/(N*K) *X11tY11 + rho11*(Gamma11 - V11)
    B11 <- solve(lhs, rhs)
    
    #### B12 update ####
    lhs = 1/(N*K) *X12tX12 + rho12*Ip
    rhs = 1/(N*K) *X12tY12 + rho12*(Gamma12 - V12)
    B12 <- solve(lhs, rhs)
    
    #### B13 update ####
    lhs = 1/(N*K) *X13tX13 + rho13*Ip
    rhs = 1/(N*K) *X13tY13 + rho13*(Gamma13 - V13)
    B13 <- solve(lhs, rhs)
    
    ####  ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)
    
    pair1 <- update_Gamma01_Gamma1_exact(B0, B1, U1, V1, lambda1, rho01, rho1, p, K)
    Gamma01 <- pair1$Gamma01
    Gamma1  <- pair1$Gamma1
    
    pair2 <- update_Gamma01_Gamma1_exact(B0, B2, U2, V2, lambda2, rho02, rho2, p, K)
    Gamma02 <- pair2$Gamma01
    Gamma2  <- pair2$Gamma1
    
    pair3 <- update_Gamma01_Gamma1_exact(B0, B3, U3, V3, lambda3, rho03, rho3, p, K)
    Gamma03 <- pair3$Gamma01
    Gamma3  <- pair3$Gamma1
    
    pair4 <- update_Gamma01_Gamma1_exact(B0, B4, U4, V4, lambda4, rho04, rho4, p, K)
    Gamma04 <- pair4$Gamma01
    Gamma4  <- pair4$Gamma1
    
    pair5 <- update_Gamma01_Gamma1_exact(B0, B5, U5, V5, lambda5, rho05, rho5, p, K)
    Gamma05 <- pair5$Gamma01
    Gamma5  <- pair5$Gamma1
    
    pair6 <- update_Gamma01_Gamma1_exact(B0, B6, U6, V6, lambda6, rho06, rho6, p, K)
    Gamma06 <- pair6$Gamma01
    Gamma6  <- pair6$Gamma1
    
    pair7 <- update_Gamma01_Gamma1_exact(B0, B7, U7, V7, lambda7, rho07, rho7, p, K)
    Gamma07 <- pair7$Gamma01
    Gamma7  <- pair7$Gamma1
    
    pair8 <- update_Gamma01_Gamma1_exact(B0, B8, U8, V8, lambda8, rho08, rho8, p, K)
    Gamma08 <- pair8$Gamma01
    Gamma8  <- pair8$Gamma1
    
    pair9 <- update_Gamma01_Gamma1_exact(B0, B9, U9, V9, lambda9, rho09, rho9, p, K)
    Gamma09 <- pair9$Gamma01
    Gamma9  <- pair9$Gamma1
    
    pair10 <- update_Gamma01_Gamma1_exact(B0, B10, U10, V10, lambda10, rho010, rho10, p, K)
    Gamma010 <- pair10$Gamma01
    Gamma10  <- pair10$Gamma1
    
    pair11 <- update_Gamma01_Gamma1_exact(B0, B11, U11, V11, lambda11, rho011, rho11, p, K)
    Gamma011 <- pair11$Gamma01
    Gamma11  <- pair11$Gamma1
    
    pair12 <- update_Gamma01_Gamma1_exact(B0, B12, U12, V12, lambda12, rho012, rho12, p, K)
    Gamma012 <- pair12$Gamma01
    Gamma12  <- pair12$Gamma1
    
    pair13 <- update_Gamma01_Gamma1_exact(B0, B13, U13, V13, lambda13, rho013, rho13, p, K)
    Gamma013 <- pair13$Gamma01
    Gamma13  <- pair13$Gamma1
    
    ####  Dual updates ####
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
    r_primal00 = Frob_norm(B0 - Gamma00) ; s_dual00 = rho00 * Frob_norm(Gamma00 - Gamma00_old)
    r_primal01 = Frob_norm(B0 - Gamma01) ; s_dual01 = rho01 * Frob_norm(Gamma01 - Gamma01_old)
    r_primal02 = Frob_norm(B0 - Gamma02) ; s_dual02 = rho02 * Frob_norm(Gamma02 - Gamma02_old)
    r_primal03 = Frob_norm(B0 - Gamma03) ; s_dual03 = rho03 * Frob_norm(Gamma03 - Gamma03_old)
    r_primal04 = Frob_norm(B0 - Gamma04) ; s_dual04 = rho04 * Frob_norm(Gamma04 - Gamma04_old)
    r_primal05 = Frob_norm(B0 - Gamma05) ; s_dual05 = rho05 * Frob_norm(Gamma05 - Gamma05_old)
    r_primal06 = Frob_norm(B0 - Gamma06) ; s_dual06 = rho06 * Frob_norm(Gamma06 - Gamma06_old)
    r_primal07 = Frob_norm(B0 - Gamma07) ; s_dual07 = rho07 * Frob_norm(Gamma07 - Gamma07_old)
    r_primal08 = Frob_norm(B0 - Gamma08) ; s_dual08 = rho08 * Frob_norm(Gamma08 - Gamma08_old)
    r_primal09 = Frob_norm(B0 - Gamma09) ; s_dual09 = rho09 * Frob_norm(Gamma09 - Gamma09_old)
    r_primal010 = Frob_norm(B0 - Gamma010) ; s_dual010 = rho010 * Frob_norm(Gamma010 - Gamma010_old)
    r_primal011 = Frob_norm(B0 - Gamma011) ; s_dual011 = rho011 * Frob_norm(Gamma011 - Gamma011_old)
    r_primal012 = Frob_norm(B0 - Gamma012) ; s_dual012 = rho012 * Frob_norm(Gamma012 - Gamma012_old)
    r_primal013 = Frob_norm(B0 - Gamma013) ; s_dual013 = rho013 * Frob_norm(Gamma013 - Gamma013_old)
    
    r_primal1 = Frob_norm(B1 - Gamma1) ; s_dual1 = rho1 * Frob_norm(Gamma1 - Gamma1_old)
    r_primal2 = Frob_norm(B2 - Gamma2) ; s_dual2 = rho2 * Frob_norm(Gamma2 - Gamma2_old)
    r_primal3 = Frob_norm(B3 - Gamma3) ; s_dual3 = rho3 * Frob_norm(Gamma3 - Gamma3_old)
    r_primal4 = Frob_norm(B4 - Gamma4) ; s_dual4 = rho4 * Frob_norm(Gamma4 - Gamma4_old)
    r_primal5 = Frob_norm(B5 - Gamma5) ; s_dual5 = rho5 * Frob_norm(Gamma5 - Gamma5_old)
    r_primal6 = Frob_norm(B6 - Gamma6) ; s_dual6 = rho6 * Frob_norm(Gamma6 - Gamma6_old)
    r_primal7 = Frob_norm(B7 - Gamma7) ; s_dual7 = rho7 * Frob_norm(Gamma7 - Gamma7_old)
    r_primal8 = Frob_norm(B8 - Gamma8) ; s_dual8 = rho8 * Frob_norm(Gamma8 - Gamma8_old)
    r_primal9 = Frob_norm(B9 - Gamma9) ; s_dual9 = rho9 * Frob_norm(Gamma9 - Gamma9_old)
    r_primal10 = Frob_norm(B10 - Gamma10) ; s_dual10 = rho10 * Frob_norm(Gamma10 - Gamma10_old)
    r_primal11 = Frob_norm(B11 - Gamma11) ; s_dual11 = rho11 * Frob_norm(Gamma11 - Gamma11_old)
    r_primal12 = Frob_norm(B12 - Gamma12) ; s_dual12 = rho12 * Frob_norm(Gamma12 - Gamma12_old)
    r_primal13 = Frob_norm(B13 - Gamma13) ; s_dual13 = rho13 * Frob_norm(Gamma13 - Gamma13_old)
    
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
  
  end <- Sys.time()
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6, B7 = B7, B8 = B8, B9 = B9, B10 = B10, B11 = B11, B12 = B12, B13 = B13,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06, Gamma07 = Gamma07, Gamma08 = Gamma08, Gamma09 = Gamma09, Gamma010 = Gamma010, Gamma011 = Gamma011, Gamma012 = Gamma012, Gamma013 = Gamma013,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6, Gamma7 = Gamma7, Gamma8 = Gamma8, Gamma9 = Gamma9, Gamma10 = Gamma10, Gamma11 = Gamma11, Gamma12 = Gamma12, Gamma13 = Gamma13,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6, U7 = U7, U8 = U8, U9 = U9, U10 = U10, U11 = U11, U12 = U12, U13 = U13,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6, V7 = V7, V8 = V8, V9 = V9, V10 = V10, V11 = V11, V12 = V12, V13 = V13,
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
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

admm_nine_source <- function(X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5, X6, Y6, X7, Y7, X8, Y8, X9, Y9,
                             X0tX0 = NULL, X0tY0 = NULL,
                             X1tX1 = NULL, X1tY1 = NULL,
                             X2tX2 = NULL, X2tY2 = NULL,
                             X3tX3 = NULL, X3tY3 = NULL,
                             X4tX4 = NULL, X4tY4 = NULL,
                             X5tX5 = NULL, X5tY5 = NULL,
                             X6tX6 = NULL, X6tY6 = NULL,
                             X7tX7 = NULL, X7tY7 = NULL,
                             X8tX8 = NULL, X8tY8 = NULL,
                             X9tX9 = NULL, X9tY9 = NULL,
                             lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8, lambda9,
                             rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1, rho07 = 1, rho08 = 1, rho09 = 1,
                             rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1, rho7 = 1, rho8 = 1, rho9 = 1,
                             mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                             max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F){
  start <- Sys.time()
  
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
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K)
  B1 <- matrix(0, p, K)
  B2 <- matrix(0, p, K)
  B3 <- matrix(0, p, K)
  B4 <- matrix(0, p, K)
  B5 <- matrix(0, p, K)
  B6 <- matrix(0, p, K)
  B7 <- matrix(0, p, K)
  B8 <- matrix(0, p, K)
  B9 <- matrix(0, p, K)
  
  Gamma00 <- matrix(0, p, K)
  Gamma01 <- matrix(0, p, K)
  Gamma02 <- matrix(0, p, K)
  Gamma03 <- matrix(0, p, K)
  Gamma04 <- matrix(0, p, K)
  Gamma05 <- matrix(0, p, K)
  Gamma06 <- matrix(0, p, K)
  Gamma07 <- matrix(0, p, K)
  Gamma08 <- matrix(0, p, K)
  Gamma09 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K)
  Gamma2 <- matrix(0, p, K)
  Gamma3 <- matrix(0, p, K)
  Gamma4 <- matrix(0, p, K)
  Gamma5 <- matrix(0, p, K)
  Gamma6 <- matrix(0, p, K)
  Gamma7 <- matrix(0, p, K)
  Gamma8 <- matrix(0, p, K)
  Gamma9 <- matrix(0, p, K)
  
  U0 <- matrix(0, p, K)
  U1 <- matrix(0, p, K)
  U2 <- matrix(0, p, K)
  U3 <- matrix(0, p, K)
  U4 <- matrix(0, p, K)
  U5 <- matrix(0, p, K)
  U6 <- matrix(0, p, K)
  U7 <- matrix(0, p, K)
  U8 <- matrix(0, p, K)
  U9 <- matrix(0, p, K)
  
  V1 <- matrix(0, p, K)
  V2 <- matrix(0, p, K)
  V3 <- matrix(0, p, K)
  V4 <- matrix(0, p, K)
  V5 <- matrix(0, p, K)
  V6 <- matrix(0, p, K)
  V7 <- matrix(0, p, K)
  V8 <- matrix(0, p, K)
  V9 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    X0tX0 <- crossprod(X0, X0)
    X1tX1 <- crossprod(X1, X1)
    X2tX2 <- crossprod(X2, X2)
    X3tX3 <- crossprod(X3, X3)
    X4tX4 <- crossprod(X4, X4)
    X5tX5 <- crossprod(X5, X5)
    X6tX6 <- crossprod(X6, X6)
    X7tX7 <- crossprod(X7, X7)
    X8tX8 <- crossprod(X8, X8)
    X9tX9 <- crossprod(X9, X9)
    
    X0tY0 <- crossprod(X0, Y0)
    X1tY1 <- crossprod(X1, Y1)
    X2tY2 <- crossprod(X2, Y2)
    X3tY3 <- crossprod(X3, Y3)
    X4tY4 <- crossprod(X4, Y4)
    X5tY5 <- crossprod(X5, Y5)
    X6tY6 <- crossprod(X6, Y6)
    X7tY7 <- crossprod(X7, Y7)
    X8tY8 <- crossprod(X8, Y8)
    X9tY9 <- crossprod(X9, Y9)
  }
  Ip <- diag(p)
  
  for (iter in 1:max_iter){ # main iteration
    
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
    lhs <- 1/(N*K) *X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06 + rho07 + rho08 + rho09) * Ip
    rhs <- 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6) + rho07 * (Gamma07 - U7) + rho08 * (Gamma08 - U8) + rho09 * (Gamma09 - U9)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/(N*K) *X1tX1 + rho1*Ip
    rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/(N*K) *X2tX2 + rho2*Ip
    rhs = 1/(N*K) *X2tY2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/(N*K) *X3tX3 + rho3*Ip
    rhs = 1/(N*K) *X3tY3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/(N*K) *X4tX4 + rho4*Ip
    rhs = 1/(N*K) *X4tY4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### B5 update ####
    lhs = 1/(N*K) *X5tX5 + rho5*Ip
    rhs = 1/(N*K) *X5tY5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    #### B6 update ####
    lhs = 1/(N*K) *X6tX6 + rho6*Ip
    rhs = 1/(N*K) *X6tY6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    
    #### B7 update ####
    lhs = 1/(N*K) *X7tX7 + rho7*Ip
    rhs = 1/(N*K) *X7tY7 + rho7*(Gamma7 - V7)
    B7 <- solve(lhs, rhs)
    
    #### B8 update ####
    lhs = 1/(N*K) *X8tX8 + rho8*Ip
    rhs = 1/(N*K) *X8tY8 + rho8*(Gamma8 - V8)
    B8 <- solve(lhs, rhs)
    
    #### B9 update ####
    lhs = 1/(N*K) *X9tX9 + rho9*Ip
    rhs = 1/(N*K) *X9tY9 + rho9*(Gamma9 - V9)
    B9 <- solve(lhs, rhs)
    
    ####  ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)
    
    pair1 <- update_Gamma01_Gamma1_exact(B0, B1, U1, V1, lambda1, rho01, rho1, p, K)
    Gamma01 <- pair1$Gamma01
    Gamma1  <- pair1$Gamma1
    
    pair2 <- update_Gamma01_Gamma1_exact(B0, B2, U2, V2, lambda2, rho02, rho2, p, K)
    Gamma02 <- pair2$Gamma01
    Gamma2  <- pair2$Gamma1
    
    pair3 <- update_Gamma01_Gamma1_exact(B0, B3, U3, V3, lambda3, rho03, rho3, p, K)
    Gamma03 <- pair3$Gamma01
    Gamma3  <- pair3$Gamma1
    
    pair4 <- update_Gamma01_Gamma1_exact(B0, B4, U4, V4, lambda4, rho04, rho4, p, K)
    Gamma04 <- pair4$Gamma01
    Gamma4  <- pair4$Gamma1
    
    pair5 <- update_Gamma01_Gamma1_exact(B0, B5, U5, V5, lambda5, rho05, rho5, p, K)
    Gamma05 <- pair5$Gamma01
    Gamma5  <- pair5$Gamma1
    
    pair6 <- update_Gamma01_Gamma1_exact(B0, B6, U6, V6, lambda6, rho06, rho6, p, K)
    Gamma06 <- pair6$Gamma01
    Gamma6  <- pair6$Gamma1
    
    pair7 <- update_Gamma01_Gamma1_exact(B0, B7, U7, V7, lambda7, rho07, rho7, p, K)
    Gamma07 <- pair7$Gamma01
    Gamma7  <- pair7$Gamma1
    
    pair8 <- update_Gamma01_Gamma1_exact(B0, B8, U8, V8, lambda8, rho08, rho8, p, K)
    Gamma08 <- pair8$Gamma01
    Gamma8  <- pair8$Gamma1
    
    pair9 <- update_Gamma01_Gamma1_exact(B0, B9, U9, V9, lambda9, rho09, rho9, p, K)
    Gamma09 <- pair9$Gamma01
    Gamma9  <- pair9$Gamma1
    
    ####  Dual updates ####
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
    r_primal00 = Frob_norm(B0 - Gamma00) ; s_dual00 = rho00 * Frob_norm(Gamma00 - Gamma00_old)
    r_primal01 = Frob_norm(B0 - Gamma01) ; s_dual01 = rho01 * Frob_norm(Gamma01 - Gamma01_old)
    r_primal02 = Frob_norm(B0 - Gamma02) ; s_dual02 = rho02 * Frob_norm(Gamma02 - Gamma02_old)
    r_primal03 = Frob_norm(B0 - Gamma03) ; s_dual03 = rho03 * Frob_norm(Gamma03 - Gamma03_old)
    r_primal04 = Frob_norm(B0 - Gamma04) ; s_dual04 = rho04 * Frob_norm(Gamma04 - Gamma04_old)
    r_primal05 = Frob_norm(B0 - Gamma05) ; s_dual05 = rho05 * Frob_norm(Gamma05 - Gamma05_old)
    r_primal06 = Frob_norm(B0 - Gamma06) ; s_dual06 = rho06 * Frob_norm(Gamma06 - Gamma06_old)
    r_primal07 = Frob_norm(B0 - Gamma07) ; s_dual07 = rho07 * Frob_norm(Gamma07 - Gamma07_old)
    r_primal08 = Frob_norm(B0 - Gamma08) ; s_dual08 = rho08 * Frob_norm(Gamma08 - Gamma08_old)
    r_primal09 = Frob_norm(B0 - Gamma09) ; s_dual09 = rho09 * Frob_norm(Gamma09 - Gamma09_old)
    
    r_primal1 = Frob_norm(B1 - Gamma1) ; s_dual1 = rho1 * Frob_norm(Gamma1 - Gamma1_old)
    r_primal2 = Frob_norm(B2 - Gamma2) ; s_dual2 = rho2 * Frob_norm(Gamma2 - Gamma2_old)
    r_primal3 = Frob_norm(B3 - Gamma3) ; s_dual3 = rho3 * Frob_norm(Gamma3 - Gamma3_old)
    r_primal4 = Frob_norm(B4 - Gamma4) ; s_dual4 = rho4 * Frob_norm(Gamma4 - Gamma4_old)
    r_primal5 = Frob_norm(B5 - Gamma5) ; s_dual5 = rho5 * Frob_norm(Gamma5 - Gamma5_old)
    r_primal6 = Frob_norm(B6 - Gamma6) ; s_dual6 = rho6 * Frob_norm(Gamma6 - Gamma6_old)
    r_primal7 = Frob_norm(B7 - Gamma7) ; s_dual7 = rho7 * Frob_norm(Gamma7 - Gamma7_old)
    r_primal8 = Frob_norm(B8 - Gamma8) ; s_dual8 = rho8 * Frob_norm(Gamma8 - Gamma8_old)
    r_primal9 = Frob_norm(B9 - Gamma9) ; s_dual9 = rho9 * Frob_norm(Gamma9 - Gamma9_old)
    
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
  
  end <- Sys.time()
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6, B7 = B7, B8 = B8, B9 = B9,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06, Gamma07 = Gamma07, Gamma08 = Gamma08, Gamma09 = Gamma09,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6, Gamma7 = Gamma7, Gamma8 = Gamma8, Gamma9 = Gamma9,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6, U7 = U7, U8 = U8, U9 = U9,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6, V7 = V7, V8 = V8, V9 = V9,
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
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

admm_six_source <- function(X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5, X6, Y6,
                            X0tX0 = NULL, X0tY0 = NULL,
                            X1tX1 = NULL, X1tY1 = NULL,
                            X2tX2 = NULL, X2tY2 = NULL,
                            X3tX3 = NULL, X3tY3 = NULL,
                            X4tX4 = NULL, X4tY4 = NULL,
                            X5tX5 = NULL, X5tY5 = NULL,
                            X6tX6 = NULL, X6tY6 = NULL,
                            lambda0, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6,
                            rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho04 = 1, rho05 = 1, rho06 = 1,
                            rho1 = 1, rho2 = 1, rho3 = 1, rho4 = 1, rho5 = 1, rho6 = 1,
                            mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                            max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F){
  start <- Sys.time()
  
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
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K)
  B1 <- matrix(0, p, K)
  B2 <- matrix(0, p, K)
  B3 <- matrix(0, p, K)
  B4 <- matrix(0, p, K)
  B5 <- matrix(0, p, K)
  B6 <- matrix(0, p, K)
  
  Gamma00 <- matrix(0, p, K)
  Gamma01 <- matrix(0, p, K)
  Gamma02 <- matrix(0, p, K)
  Gamma03 <- matrix(0, p, K)
  Gamma04 <- matrix(0, p, K)
  Gamma05 <- matrix(0, p, K)
  Gamma06 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K)
  Gamma2 <- matrix(0, p, K)
  Gamma3 <- matrix(0, p, K)
  Gamma4 <- matrix(0, p, K)
  Gamma5 <- matrix(0, p, K)
  Gamma6 <- matrix(0, p, K)
  
  U0 <- matrix(0, p, K)
  U1 <- matrix(0, p, K)
  U2 <- matrix(0, p, K)
  U3 <- matrix(0, p, K)
  U4 <- matrix(0, p, K)
  U5 <- matrix(0, p, K)
  U6 <- matrix(0, p, K)
  
  V1 <- matrix(0, p, K)
  V2 <- matrix(0, p, K)
  V3 <- matrix(0, p, K)
  V4 <- matrix(0, p, K)
  V5 <- matrix(0, p, K)
  V6 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    X0tX0 <- crossprod(X0, X0)
    X1tX1 <- crossprod(X1, X1)
    X2tX2 <- crossprod(X2, X2)
    X3tX3 <- crossprod(X3, X3)
    X4tX4 <- crossprod(X4, X4)
    X5tX5 <- crossprod(X5, X5)
    X6tX6 <- crossprod(X6, X6)
    
    X0tY0 <- crossprod(X0, Y0)
    X1tY1 <- crossprod(X1, Y1)
    X2tY2 <- crossprod(X2, Y2)
    X3tY3 <- crossprod(X3, Y3)
    X4tY4 <- crossprod(X4, Y4)
    X5tY5 <- crossprod(X5, Y5)
    X6tY6 <- crossprod(X6, Y6)
  }
  Ip <- diag(p)
  
  for (iter in 1:max_iter){ # main iteration
    
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
    lhs <- 1/(N*K) *X0tX0 + (rho00 + rho01 + rho02 + rho03 + rho04 + rho05 + rho06) * Ip
    rhs <- 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3) + rho04 * (Gamma04 - U4) + rho05 * (Gamma05 - U5) + rho06 * (Gamma06 - U6)
    B0 <- solve(lhs, rhs)
    
    #### B1 update ####
    lhs = 1/(N*K) *X1tX1 + rho1*Ip
    rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    B1 <- solve(lhs, rhs)
    
    #### B2 update ####
    lhs = 1/(N*K) *X2tX2 + rho2*Ip
    rhs = 1/(N*K) *X2tY2 + rho2*(Gamma2 - V2)
    B2 <- solve(lhs, rhs)
    
    #### B3 update ####
    lhs = 1/(N*K) *X3tX3 + rho3*Ip
    rhs = 1/(N*K) *X3tY3 + rho3*(Gamma3 - V3)
    B3 <- solve(lhs, rhs)
    
    #### B4 update ####
    lhs = 1/(N*K) *X4tX4 + rho4*Ip
    rhs = 1/(N*K) *X4tY4 + rho4*(Gamma4 - V4)
    B4 <- solve(lhs, rhs)
    
    #### B5 update ####
    lhs = 1/(N*K) *X5tX5 + rho5*Ip
    rhs = 1/(N*K) *X5tY5 + rho5*(Gamma5 - V5)
    B5 <- solve(lhs, rhs)
    
    #### B6 update ####
    lhs = 1/(N*K) *X6tX6 + rho6*Ip
    rhs = 1/(N*K) *X6tY6 + rho6*(Gamma6 - V6)
    B6 <- solve(lhs, rhs)
    
    ####  ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)
    
    pair1 <- update_Gamma01_Gamma1_exact(B0, B1, U1, V1, lambda1, rho01, rho1, p, K)
    Gamma01 <- pair1$Gamma01
    Gamma1  <- pair1$Gamma1
    
    pair2 <- update_Gamma01_Gamma1_exact(B0, B2, U2, V2, lambda2, rho02, rho2, p, K)
    Gamma02 <- pair2$Gamma01
    Gamma2  <- pair2$Gamma1
    
    pair3 <- update_Gamma01_Gamma1_exact(B0, B3, U3, V3, lambda3, rho03, rho3, p, K)
    Gamma03 <- pair3$Gamma01
    Gamma3  <- pair3$Gamma1
    
    pair4 <- update_Gamma01_Gamma1_exact(B0, B4, U4, V4, lambda4, rho04, rho4, p, K)
    Gamma04 <- pair4$Gamma01
    Gamma4  <- pair4$Gamma1
    
    pair5 <- update_Gamma01_Gamma1_exact(B0, B5, U5, V5, lambda5, rho05, rho5, p, K)
    Gamma05 <- pair5$Gamma01
    Gamma5  <- pair5$Gamma1
    
    pair6 <- update_Gamma01_Gamma1_exact(B0, B6, U6, V6, lambda6, rho06, rho6, p, K)
    Gamma06 <- pair6$Gamma01
    Gamma6  <- pair6$Gamma1
    
    ####  Dual updates ####
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
    r_primal00 = Frob_norm(B0 - Gamma00) ; s_dual00 = rho00 * Frob_norm(Gamma00 - Gamma00_old)
    r_primal01 = Frob_norm(B0 - Gamma01) ; s_dual01 = rho01 * Frob_norm(Gamma01 - Gamma01_old)
    r_primal02 = Frob_norm(B0 - Gamma02) ; s_dual02 = rho02 * Frob_norm(Gamma02 - Gamma02_old)
    r_primal03 = Frob_norm(B0 - Gamma03) ; s_dual03 = rho03 * Frob_norm(Gamma03 - Gamma03_old)
    r_primal04 = Frob_norm(B0 - Gamma04) ; s_dual04 = rho04 * Frob_norm(Gamma04 - Gamma04_old)
    r_primal05 = Frob_norm(B0 - Gamma05) ; s_dual05 = rho05 * Frob_norm(Gamma05 - Gamma05_old)
    r_primal06 = Frob_norm(B0 - Gamma06) ; s_dual06 = rho06 * Frob_norm(Gamma06 - Gamma06_old)
    
    r_primal1 = Frob_norm(B1 - Gamma1) ; s_dual1 = rho1 * Frob_norm(Gamma1 - Gamma1_old)
    r_primal2 = Frob_norm(B2 - Gamma2) ; s_dual2 = rho2 * Frob_norm(Gamma2 - Gamma2_old)
    r_primal3 = Frob_norm(B3 - Gamma3) ; s_dual3 = rho3 * Frob_norm(Gamma3 - Gamma3_old)
    r_primal4 = Frob_norm(B4 - Gamma4) ; s_dual4 = rho4 * Frob_norm(Gamma4 - Gamma4_old)
    r_primal5 = Frob_norm(B5 - Gamma5) ; s_dual5 = rho5 * Frob_norm(Gamma5 - Gamma5_old)
    r_primal6 = Frob_norm(B6 - Gamma6) ; s_dual6 = rho6 * Frob_norm(Gamma6 - Gamma6_old)
    
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
  
  end <- Sys.time()
  
  return(list(B0 = B0, B1 = B1, B2 = B2, B3 = B3, B4 = B4, B5 = B5, B6 = B6,
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma04 = Gamma04, Gamma05 = Gamma05, Gamma06 = Gamma06,
              Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3, Gamma4 = Gamma4, Gamma5 = Gamma5, Gamma6 = Gamma6,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, U4 = U4, U5 = U5, U6 = U6,
              V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6,
              iter = iter, runtime = as.numeric(difftime(end, start, units = "secs")),
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
