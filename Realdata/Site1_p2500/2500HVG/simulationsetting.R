# lambda_values = 10^seq(-4, 0, length = 7)
# 
# train_prop = 0.5 ; val_prop = 0.25 ; test_prop = 0.25 ; nrep = 2
# 
# tol_vanilla = 1e-3
# max_iter_vanilla = 50
# tol_prim = 1e-3
# tol_dual = 1e-3
# max_iter = 200


lambda_values = 10^seq(-6, 0, length = 15)

train_prop = 0.5 ; val_prop = 0.25 ; test_prop = 0.25 ; nrep = 20

K = 77

max_iter_vanilla_merged = 1000
tol_vanilla_merged = sqrt(K)*1e-4


max_iter_vanilla = 1000  
tol_vanilla = sqrt(K)*1e-4


max_iter = 1000
tol_prim = sqrt(K)*1e-4
tol_dual = sqrt(K)*1e-4

max_iter_stl = 200
tol_prim_stl = 5e-3
tol_dual_stl = 5e-3