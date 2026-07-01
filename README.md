# TMTL

Transfer Multitask Learning for high-dimensional multi-output regression.

## Quick start
- Generate example data
- Fit TMTL(Fused)
- Fit TMTL(Debiased)

## Sample data generation (with 2 source domains; $L=2$)

Example data generated in the same way as "Balanced shift" in the manuscript.

```r
set.seed(1)
L=2
p = 100 ; K = 50 ; n_source = 500 ; n_target = 100 ; N = L*n_source + n_target
lambda_values = 10^seq(-4, -1, length = 5)
cv_folds = 3

active_features_target = 1:10 ; active_features_source = 1:15
sigma = 0.25 ; sddiff = 0.1 ; sd_target = 1

proportion_identical_tasks = 0
Wmean = 0.2

W_true_target = matrix(0, p, K) ; Delta_W = matrix(0, p, K)
for(i in 1:length(active_features_target)){
  W_true_target[i,] = rnorm(K, mean = Wmean, sd = Wmean) # /sqrt(i), sd = Wmean/sqrt(i))
  print(norm(W_true_target[i,],"2"))
}
alpha = 2

for(i in 1:length(active_features_source)){
  Delta_W[i,] = rnorm(K, mean = Wmean/(alpha*sqrt(i)), sd = Wmean/(alpha*sqrt(i))) 
  print(norm(Delta_W[i,],"2"))
}
if(proportion_identical_tasks > 0){
  Delta_W[,(K*(1-proportion_identical_tasks)+1):K] = 0
}

W_true_source1 = W_true_target + 1*Delta_W
W_true_source2 = W_true_target - 1*Delta_W
```


## TMTL(Fused)

The following is an example code to fit \textbf{TMTL(Fused)} with $\lambda_0 = \lambda_1 = \lambda_2 = 0.001$.


```r
lambda0 = 0.001
lambda1 = 0.001
lambda2 = 0.001

set.seed(1) # seed to generate errors


## generate target and source data
X_target = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) # ; X_target = scale(X_target, center = T, scale = F) #; X_target_center = attr(X_target, "scaled:center")
  X_target_test = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) # iid copy of X_target to evaluate out of sample RMSE.
  Y_target = X_target %*% W_true_target + matrix(rnorm(n_target * K, mean = 0, sd = sigma), n_target, K) #; Y_target = scale(Y_target, center = T, scale = F) ; Y_target_center = attr(Y_target, "scaled:center")

X_source1 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) # ; X_source1 = scale(X_source1, center = T, scale = F) ; X_source1_center = attr(X_source1, "scaled:center")
  Y_source1 = X_source1 %*% W_true_source1 + matrix(rnorm(n_source * K, sd = sigma), n_source, K) # ; Y_source1 = scale(Y_source1, center = T, scale = F) ; Y_source1_center = attr(Y_source1, "scaled:center")
  X_source2 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) # ; X_source2 = scale(X_source2, center = T, scale = F) ; X_source2_center = attr(X_source2, "scaled:center")
  Y_source2 = X_source2 %*% W_true_source2 + matrix(rnorm(n_source * K, sd = sigma), n_source, K) # ; Y_source2 = scale(Y_source2, center = T, scale = F) ; Y_source2_center = attr(Y_source2, "scaled:center")
```



### precompute matrices
X0tX0 = crossprod(X_target,X_target) ; X0tY0 = crossprod(X_target,Y_target)
X1tX1 = crossprod(X_source1,X_source1) ; X1tY1 = crossprod(X_source1,Y_source1)
X2tX2 = crossprod(X_source2,X_source2) ; X2tY2 = crossprod(X_source2,Y_source2)
  

### run algorithm
admmm = admm_two_source(X0 = X_target, Y0 = Y_target, X1 = X_source1, Y1 = Y_source1, X2 = X_source2, Y2 = Y_source2,
X0tX0, X0tY0, X1tX1, X1tY1, X2tX2, X2tY2,
lambda0 = lambda0, lambda1 = lambda1, lambda2 = lambda2,
rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1,
mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)

W_opt_fused = n_source/N*(admmm$Gamma1 + admmm$Gamma2) + n_target/N*admmm$Gamma00
```



