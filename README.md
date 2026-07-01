# TMTL

Transfer Multitask Learning for high-dimensional multi-output regression.

## Table of Contents
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

First, solve the jointly regularized least-squares problem
$\hat{\mathbf{B}}^0,\cdots,\hat{\mathbf{B}}^L \in \arg\min_{\mathbf{B}^0,\cdots,\mathbf{B}^L} \left[
\frac{1}{2NK}\sum\limits_{\ell=0}^{L} \lVert \mathbf{Y}^{\ell} - \mathbf{X}^{\ell \cdot} \mathbf{B}^\ell \rVert_F^2 + \lambda_0
\left( \lVert \mathbf{B}^{0} \rVert_{2,1} + \sum\limits_{\ell=1}^{L} a_\ell \lVert \mathbf{B}^\ell - \mathbf{B}^0 \rVert_{2,1} \right) \right]$

Then compute the fused estimator:
$\mathbf{W}^F = \frac{n_S}{N} \sum_{\ell=1}^{L} \hat{\mathbf{B}}^\ell + \frac{n_T}{N} \hat{\mathbf{B}}^0$
The following is an example code to fit **TMTL(Fused)** with $\lambda_0 = \lambda_1 = \lambda_2 = 0.001$.

First, set tuning parameters and generate data.

```r
lambda0 = 0.001
lambda1 = 0.001
lambda2 = 0.001

set.seed(1) # seed to generate errors


## generate target and source data
X_target = matrix(rnorm(n_target * p, mean = 0, sd = sd_target), n_target, p) 
Y_target = X_target %*% W_true_target + matrix(rnorm(n_target * K, mean = 0, sd = sigma), n_target, K) 

X_source1 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) 
Y_source1 = X_source1 %*% W_true_source1 + matrix(rnorm(n_source * K, sd = sigma), n_source, K) 
X_source2 = matrix(rnorm(n_source * p, mean = 0, sd = sample(c(sd_target + sddiff, sd_target - sddiff), size = n_source * p, replace = T, prob = c(0.5, 0.5))), n_source, p) 
Y_source2 = X_source2 %*% W_true_source2 + matrix(rnorm(n_source * K, sd = sigma), n_source, K) 
```

Precomputing gram matrices and saving those help reduce computational cost.

```r
X0tX0 = crossprod(X_target,X_target) ; X0tY0 = crossprod(X_target,Y_target)
X1tX1 = crossprod(X_source1,X_source1) ; X1tY1 = crossprod(X_source1,Y_source1)
X2tX2 = crossprod(X_source2,X_source2) ; X2tY2 = crossprod(X_source2,Y_source2)
```

Now run the algorithm with alternating direction method of multipliers (ADMM), with adaptive update rule. 
Maximum iteration = 1000 with tolerance level for dual parameters and primal parameters were set as 1e-5.
$\mathbf{W}^F = \frac{n_S}{N} \sum_{\ell=1}^{L} \hat{\mathbf{B}}^\ell + \frac{n_T}{N} \hat{\mathbf{B}}^0$.

```r
admmm = admm_two_source(X0 = X_target, Y0 = Y_target, X1 = X_source1, Y1 = Y_source1, X2 = X_source2, Y2 = Y_source2, X0tX0, X0tY0, X1tX1, X1tY1, X2tX2, X2tY2,
                        lambda0 = lambda0, lambda1 = lambda1, lambda2 = lambda2,
                        rho00 = 1, rho01 = 1, rho02 = 1, rho1 = 1, rho2 = 1,
                         mu = 10, tau_incr = 2, tau_decr = 1.5, rho_max = 10000,
                         max_iter = 1000, tol_prim = 1e-5, tol_dual = 1e-5, verbose = F)

W_opt_fused = n_source/N*(admmm$Gamma1 + admmm$Gamma2) + n_target/N*admmm$Gamma00
```


## TMTL(Debiased)

A debiasing step only involves target-data and corrects aligned shifts from target.
$\hat{\mathbf{D}} \in \argminA_{\mathbf{D} \in \mathbb{R}^{p \times K}} \frac{1}{2n_T K} \lVert \mathbf{Y}^0 - \mathbf{X}^{0 \cdot} (\hat{\mathbf{W}}^{F} + \mathbf{D}) \rVert_F^2 + \tilde{\lambda} \lVert \mathbf{D} \rVert_{2,1}, $

and the final debiased estimator is
$\hat{\mathbf{W}}^{D} := \hat{\mathbf{W}}^{F} + \hat{\mathbf{D}}.$

The following is the code to fit **TMTL(Debiased)** with $\tilde{\lambda} = 0.001$.




