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



