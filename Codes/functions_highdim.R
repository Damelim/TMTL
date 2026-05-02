library(Rcpp)
library(RcppArmadillo)
library(Matrix) 
#sourceCpp('~/Library/CloudStorage/GoogleDrive-96limtotoro@gmail.com/My Drive/research_multitask/multitask_supplement.cpp')
update_Gamma00 <- function(B0, U0, lambda0, rho00,p,K) {
  Z <- B0 + U0
  row_norms <- sqrt(rowSums(Z^2))
  shrink <- pmax(0, 1 - lambda0 / (rho00 * row_norms))
  shrink[!is.finite(shrink)] <- 0
  Gamma00 <- Z * matrix(shrink, nrow = p, ncol = K)
  return(Gamma00)
}

update_Gamma01 <- function(B0, U1, Gamma1, lambda1, rho01,p,K) {
  Z <- B0 + U1
  Diff <- Z - Gamma1
  row_norms <- sqrt(rowSums(Diff^2))
  shrink <- pmax(0, 1 - lambda1 / (rho01 * row_norms))
  shrink[!is.finite(shrink)] <- 0
  Gamma01 <- Gamma1 + Diff * matrix(shrink, nrow = p, ncol = K)
  return(Gamma01)
}

update_Gamma1 <- function(B1, V1, Gamma01, lambda1, rho1,p,K) {
  Z <- B1 + V1
  Diff <- Z - Gamma01
  row_norms <- sqrt(rowSums(Diff^2))
  shrink <- pmax(0, 1 - lambda1 / (rho1 * row_norms))
  shrink[!is.finite(shrink)] <- 0
  Gamma1 <- Gamma01 + Diff * matrix(shrink, nrow = p, ncol = K)
  return(Gamma1)
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
                    mu = 10, tau_inc = 2, tau_dec = 1.5, sparse = F,
                    verbose = FALSE) {
  
  objectivevalue = rep(NA, max_iter)
  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(Y)
  
  if(sparse){
    X = as(X, "dgCMatrix")
  }
  Xt = t(X)
  
  B <- matrix(0, p, K)
  Gamma <- matrix(0, p, K)
  U <- matrix(0, p, K)  # scaled dual variable
  
  XtY <- crossprod(X, Y)
  if(p<n){
    XtX <- crossprod(X, X)
    Ip <- Diagonal(p)
  }else{
    XXt <- X %*% Xt
    In <- Diagonal(n)
  }
  
  for (iter in 1:max_iter) {
    
    #if(verbose){if(iter %% 10 == 0){print(iter)}}

    Gamma_old <- Gamma
    
    # B-update (ridge-like step)
    if(p < n){
      A = XtX/(n*K) + rho*Ip
      cf = Cholesky(A, LDL = F)
      B = solve(cf, XtY/(n*K) + rho*(Gamma-U))
    }else{ # use woodbury matrix identity
      A = (n*K)*In + (1/rho)*XXt
      cf = Cholesky(A, LDL=F)
      RHS = XtY/(n*K) + rho * (Gamma - U)
      B = (RHS / rho) - (Xt %*% solve(cf, X %*% (RHS/(rho^2))))
    }
    #B <- solve(1/(n*K) * XtX + rho * Ip, 1/(n*K) * XtY + rho * (Gamma - U))
    
    
    
    # Gamma-update (prox for \ell_2 norm: https://math.stackexchange.com/questions/2190885/proximal-operator-of-the-euclidean-norm-l-2-norm)
    
    Z = B+U
    row_norms = sqrt(rowSums(Z^2))
    shrink <- pmax(0, 1 - lambda/(rho*row_norms))
    shrink[!is.finite(shrink)] <- 0
    Gamma <- Z * matrix(shrink, nrow=p, ncol=K)
    # Z <- B + U
    # for (j in 1:p) {
    #   z_j <- Z[j, ]
    #   norm_z <- sqrt(sum(z_j^2))
    #   if (norm_z > lambda / rho) {
    #     Gamma[j, ] <- (1 - lambda / (rho * norm_z)) * z_j
    #   } else {
    #     Gamma[j, ] <- 0
    #   }
    # }
    
    # Dual update
    U <- U + B - Gamma
    
    # Residuals
    primal_resid <- Frob_norm(B - Gamma)
    dual_resid <- rho * Frob_norm(Gamma - Gamma_old)
    
    if (verbose && iter %% 10 == 0) {
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
    
    #objectivevalue[iter] = 1/(2*n)*Frob_norm(Y - X %*% Gamma) + lambda * twoone_norm(Gamma)
    
    # Check convergence
    if (iter > 2 & primal_resid < tol && dual_resid < tol) {
      if (verbose) cat(sprintf("Converged at iter %d\n", iter))
      break
    }
  }
  return(list(sol = Gamma))
  #return(list(objective = objectivevalue[!is.na(objectivevalue)], sol = Gamma))
}

admm_three_source <- function(X0, Y0, X1, Y1, X2, Y2, X3, Y3,
                              X0tX0 = NULL, X0tY0 = NULL, X1tX1 = NULL, X1tY1 = NULL, X2tX2 = NULL, X2tY2 = NULL, X3tX3 = NULL, X3tY3 = NULL,
                              X0X0t = NULL, X1X1t = NULL, X2X2t = NULL, X3X3t = NULL, Ip = NULL, In0 = NULL, In1 = NULL, In2 = NULL, In3 = NULL,
                              lambda0, lambda1, lambda2, lambda3,
                              rho00 = 1, rho01 = 1, rho02 = 1, rho03 = 1, rho1 = 1, rho2 = 1, rho3 = 1,
                              mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                              max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, sparse = F, verbose = F){
  
  n0 = nrow(X0) ; n1 = nrow(X1) ; n2 = nrow(X2) ; n3 = nrow(X3)
  N = n0+n1+n2+n3
  
  #objvalue = rep(NA, max_iter)
  rho00value = rep(NA, max_iter) ; rho01value = rep(NA, max_iter) ; rho02value = rep(NA, max_iter) ; rho03value = rep(NA, max_iter)
  rho1value = rep(NA, max_iter) ; rho2value = rep(NA, max_iter) ; rho3value = rep(NA, max_iter)
  
  r_primal00value = rep(NA, max_iter) ; s_dual00value = rep(NA, max_iter)
  r_primal01value = rep(NA, max_iter) ; s_dual01value = rep(NA, max_iter)
  r_primal02value = rep(NA, max_iter) ; s_dual02value = rep(NA, max_iter)
  r_primal03value = rep(NA, max_iter) ; s_dual03value = rep(NA, max_iter)
  r_primal1value = rep(NA, max_iter) ; s_dual1value = rep(NA, max_iter)
  r_primal2value = rep(NA, max_iter) ; s_dual2value = rep(NA, max_iter)
  r_primal3value = rep(NA, max_iter) ; s_dual3value = rep(NA, max_iter)
  
  if(sparse){
    X0 = as(X0, "dgCMatrix") ; X1 = as(X1, "dgCMatrix") ; X2 = as(X2, "dgCMatrix") ; X3 = as(X3, "dgCMatrix")
  }
  
  p <- ncol(X0)
  K <- ncol(Y0)
  
  B0 <- matrix(0, p, K) ; B1 <- matrix(0, p, K) ; B2 <- matrix(0, p, K) ; B3 <- matrix(0, p, K)
  Gamma00 <- matrix(0, p, K) ; Gamma01 <- matrix(0, p, K) ; Gamma02 <- matrix(0, p, K) ; Gamma03 <- matrix(0, p, K)
  Gamma1 <- matrix(0, p, K) ; Gamma2 <- matrix(0, p, K) ; Gamma3 <- matrix(0, p, K)
  U0 <- matrix(0, p, K) ; U1 <- matrix(0, p, K) ; U2 <- matrix(0, p, K) ; U3 <- matrix(0, p, K)
  V1 <- matrix(0, p, K) ; V2 <- matrix(0, p, K) ; V3 <- matrix(0, p, K)
  
  if(is.null(X0tX0)){
    if(p < min(c(n0,n1,n2,n3))){
      X0tX0 <- crossprod(X0, X0)
      X1tX1 <- crossprod(X1, X1)
      X2tX2 <- crossprod(X2, X2)
      X3tX3 <- crossprod(X3, X3)
      Ip <- Diagonal(p)
    }
  }
  if(is.null(X0tY0)){
    X0tY0 <- crossprod(X0, Y0)
    X1tY1 <- crossprod(X1, Y1)
    X2tY2 <- crossprod(X2, Y2)
    X3tY3 <- crossprod(X3, Y3)
  }
  
  X0t = t(X0) ; X1t = t(X1) ; X2t = t(X2) ; X3t = t(X3) # have to be out of the next if statement b/c X0t,... are not ftn arguments.
  if(is.null(X0X0t)){
    X0X0t = X0 %*% X0t ; X1X1t = X1 %*% X1t ; X2X2t = X2 %*% X2t ; X3X3t = X3 %*% X3t
    In0 = Diagonal(n0) ; In1 = Diagonal(n1) ; In2 = Diagonal(n2) ; In3 = Diagonal(n3)
  }

  
  
  for (iter in 1:max_iter){ # main iteration
    
    # if(verbose){
    #   if(iter %% 10 == 0){print(iter)}
    # }
    
    Gamma00_old = Gamma00 ; Gamma01_old = Gamma01 ; Gamma02_old = Gamma02 ; Gamma03_old = Gamma03
    Gamma1_old = Gamma1 ; Gamma2_old = Gamma2 ; Gamma3_old = Gamma3
    
    #### Step 1: B0 update ####
    # lhs <- 1/(N*K) *X0tX0 + (rho00 + rho01 + rho02 + rho03) * Ip
    # rhs <- 1/(N*K) *X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1) + rho02 * (Gamma02 - U2) + rho03 * (Gamma03 - U3)
    # B0 <- solve(lhs, rhs)
    
    rho0sum <- rho00 + rho01 + rho02 + rho03
    if (p < n0) {
      RHS0 <- (X0tY0 / (N*K)) + rho00*(Gamma00 - U0) + rho01*(Gamma01 - U1) + rho02*(Gamma02 - U2) + rho03*(Gamma03 - U3)
      A0  <- (X0tX0 / (N*K)) + rho0sum * Ip
      cf0 <- Cholesky(A0, LDL=FALSE)
      B0 <- solve(cf0, RHS0)
    } else {
      RHS0 <- (X0tY0 / (N*K)) + rho00*(Gamma00 - U0) + rho01*(Gamma01 - U1) + rho02*(Gamma02 - U2) + rho03*(Gamma03 - U3)
      M0  <- (N*K) * In0 + (1/rho0sum) * X0X0t
      cf0W <- Cholesky(M0, LDL=FALSE)
      B0 <- RHS0 / rho0sum - (X0t %*% solve(cf0W, X0 %*% (RHS0 / (rho0sum^2))))
    }
    
    #### Step 2: B1 update ####
    # lhs = 1/(N*K) *X1tX1 + rho1*Ip
    # rhs = 1/(N*K) *X1tY1 + rho1*(Gamma1 - V1)
    # B1 <- solve(lhs, rhs)
    
    if (p < n1) {
      RHS1 <- (X1tY1 / (N*K)) + rho1*(Gamma1 - V1)
      A1  <- (X1tX1 / (N*K)) + rho1 * Ip
      cf1 <- Cholesky(A1, LDL=FALSE)
      B1 <- solve(cf1, RHS1)
    } else {
      RHS1 <- (X1tY1 / (N*K)) + rho1*(Gamma1 - V1)
      M1  <- (N*K) * In1 + (1/rho1) * X1X1t
      cf1W <- Cholesky(M1, LDL=FALSE)
      B1 <- RHS1 / rho1 - (X1t %*% solve(cf1W, X1 %*% (RHS1 / (rho1^2))))
    }
    
    
    
    #### Step 3: B2 update ####
    # lhs = 1/(N*K) *X2tX2 + rho2*Ip
    # rhs = 1/(N*K) *X2tY2 + rho2*(Gamma2 - V2)
    # B2 <- solve(lhs, rhs)
    
    if (p < n2) {
      RHS2 <- (X2tY2 / (N*K)) + rho2*(Gamma2 - V2)
      A2  <- (X2tX2 / (N*K)) + rho2 * Ip
      cf2 <- Cholesky(A2, LDL=FALSE)
      B2 <- solve(cf2, RHS2)
    } else {
      RHS2 <- (X2tY2 / (N*K)) + rho2*(Gamma2 - V2)
      M2  <- (N*K) * In2 + (1/rho2) * X2X2t
      cf2W <- Cholesky(M2, LDL=FALSE)
      B2 <- RHS2 / rho2 - (X2t %*% solve(cf2W, X2 %*% (RHS2 / (rho2^2))))
    }
    
    
    #### Step 4: B3 update ####
    # lhs = 1/(N*K) *X3tX3 + rho3*Ip
    # rhs = 1/(N*K) *X3tY3 + rho3*(Gamma3 - V3)
    # B3 <- solve(lhs, rhs)
    
    if (p < n3) {
      RHS3 <- (X3tY3 / (N*K)) + rho3*(Gamma3 - V3)
      A3  <- (X3tX3 / (N*K)) + rho3 * Ip
      cf3 <- Cholesky(A3, LDL=FALSE)
      B3 <- solve(cf3, RHS3)
    } else {
      RHS3 <- (X3tY3 / (N*K)) + rho3*(Gamma3 - V3)
      M3  <- (N*K) * In3 + (1/rho3) * X3X3t
      cf3W <- Cholesky(M3, LDL=FALSE)
      B3 <- RHS3 / rho3 - (X3t %*% solve(cf3W, X3 %*% (RHS3 / (rho3^2))))
    }
    
    #### Step 5: Gamma00 update ####
    Gamma00 = update_Gamma00(B0,U0,lambda0,rho00,p,K)
    #### Step 6-8: Gamma01 - Gamma03 update ####
    Gamma01 = update_Gamma01(B0,U1,Gamma1,lambda1,rho01,p,K)
    Gamma02 = update_Gamma01(B0,U2,Gamma2,lambda2,rho02,p,K)
    Gamma03 = update_Gamma01(B0,U3,Gamma3,lambda3,rho03,p,K)
    
    #### Step 9-11: Gamma1 update ####
    Gamma1 = update_Gamma1(B1,V1,Gamma01,lambda1,rho1,p,K)
    Gamma2 = update_Gamma1(B2,V2,Gamma02,lambda2,rho2,p,K)
    Gamma3 = update_Gamma1(B3,V3,Gamma03,lambda3,rho3,p,K)
    
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
    
    if (verbose && iter %% 10 == 0) {
      cat(sprintf("Iter %d | Primal Resid00: %.1f | Dual Resid00: %.1f | Primal Resid01: %.1f | Dual Resid01: %.1f | Primal Resid02: %.1f | Dual Resid02: %.1f | Primal Resid03: %.1f | Dual Resid03: %.1f | Primal Resid1: %.1f | Dual Resid1: %.1f | Primal Resid2: %.1f | Dual Resid2: %.1f | Primal Resid3: %.1f | Dual Resid3: %.1f | \n",
                  iter, log(r_primal00,base=10), log(s_dual00,base=10), log(r_primal01,base=10), log(s_dual01,base=10), log(r_primal02,base=10), log(s_dual02,base=10), log(r_primal03,base=10), log(s_dual03,base=10), log(r_primal1,base=10), log(s_dual1,base=10), log(r_primal2,base=10), log(s_dual2,base=10), log(r_primal3,base=10), log(s_dual3,base=10)))
    }
    
    #objvalue[iter] = 1/(2*N) * (Frob_norm_sq(Y_target - X_target %*% B0) + Frob_norm_sq(Y_source1 - X_source1 %*% B1) + Frob_norm_sq(Y_source2 - X_source2 %*% B2) + Frob_norm_sq(Y_source3 - X_source3 %*% B3)) + lambda0 * twoone_norm(B0) + lambda1 * twoone_norm(B0 - B1) + lambda2 * twoone_norm(B0 - B2) + lambda3 * twoone_norm(B0 - B3) # objective value for ADMM
    
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
              Gamma00 = Gamma00, Gamma01 = Gamma01, Gamma02 = Gamma02, Gamma03 = Gamma03, Gamma1 = Gamma1, Gamma2 = Gamma2, Gamma3 = Gamma3,
              U0 = U0, U1 = U1, U2 = U2, U3 = U3, V1 = V1, V2 = V2, V3 = V3,
              #objvalue = objvalue, 
              rho00 = rho00value, rho01 = rho01value, rho02 = rho02value, rho03 = rho03value, rho1 = rho1value, rho2 = rho2value, rho3 = rho3value,
              r_primal00 = r_primal00value, s_dual00 = s_dual00value, 
              r_primal01 = r_primal01value, s_dual01 = s_dual01value,
              r_primal02 = r_primal02value, s_dual02 = s_dual02value,
              r_primal03 = r_primal03value, s_dual03 = s_dual03value,
              r_primal1 = r_primal1value, s_dual1 = s_dual1value,
              r_primal2 = r_primal2value, s_dual2 = s_dual2value,
              r_primal3 = r_primal3value, s_dual3 = s_dual3value))
}
