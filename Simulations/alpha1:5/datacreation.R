nreplications = 10

set.seed(1)
L=1
p = 100 ; K = 50 ; n_source = 500 ; n_target = 100 ; N = L*n_source + n_target
lambda_values = 10^seq(-4, 0, length = 7)
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
alpha = 5

for(i in 1:length(active_features_source)){
  Delta_W[i,] = rnorm(K, mean = Wmean/(alpha*sqrt(i)), sd = Wmean/(alpha*sqrt(i))) 
  print(norm(Delta_W[i,],"2"))
}
if(proportion_identical_tasks > 0){
  Delta_W[,(K*(1-proportion_identical_tasks)+1):K] = 0
}

W_true_source = W_true_target + Delta_W

print(round(W_true_target,3))
print(round(Delta_W,3))

save.image('onesource.Rdata')
rm(list=ls())