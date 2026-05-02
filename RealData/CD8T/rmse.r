##### Nofiltered 1000HVG for 5 replications X log1p --> scaled, Y CLR and not scaled afterwards (Donor "15078" only, sites 1 to 4) ####
colvec = c('red1','pink2',
           #'orange1','orange4',
           'green1','green4'
           ,
           'blue','skyblue3'
)
fname_vec = c('MTL_transfer','MTL_transfer_debiased',
              #'MTL_adaptive_transfer','MTL_adaptive_transfer_debiased',
              'MTL_notransfer','MTL_notransfer_merged'
              ,'STL_transfer_cluster','STL_transfer_debiased'
)
legendvec = c('TMTL(Fused)','TMTL(Debiased)',
              #'TMTL(Fused-adaptive)','TMTL(Debiased-adaptive)',
              'MTL(Target)','MTL(Full)'
              ,'TSTL(Fused)','TSTL(Debiased)'
)
donorvec = 1:8 
n_donor_vec = rep(NA,length(donorvec))
load('clr_data_donor_nofiltered_1000hvg1.Rdata') ; K = ncol(Y) 
for(jj in 1:length(donorvec)){
  load(paste('clr_data_donor_nofiltered_1000hvg',jj,'.Rdata', sep = ""))
  n_donor_vec[jj] = nrow(Y)
}
n_donor_vec

mse_mat_list = list()
load(paste(fname_vec[3],'/summary_logtransform_target',donorvec[1],'.Rdata',sep=""))
nrep = length(prediction_mse_vec) ; print(nrep)
for(l in 1:length(n_donor_vec)){
  cat('donor ',l,'\n')
  mse_mat_l = matrix(NA, length(fname_vec),nrep)
  for(a in 1:length(fname_vec)){
    load(paste(fname_vec[a],'/summary_logtransform_target',donorvec[l],'.Rdata',sep=""))
    if(a != 5){
      mse_mat_l[a,1:nrep] = prediction_mse_vec[1:nrep]/K
    }else if (a == 5){
      mse_mat_l[a,1:nrep] = prediction_mse_vec[1:nrep]^2 * n_donor_vec[l]/4 # divide by 4 b/c of the test data proportion 25%
    }
  }
  mse_mat_list[[l]] = mse_mat_l
}
mse_mat_list[[1]]

# for(j in 1:length(n_donor_vec)){print(rowMeans(mse_mat_list[[j]])[3] - rowMeans(mse_mat_list[[j]])[2]) # positive --> debiased better than notransfer }

mse_mat_overall = matrix(0, length(fname_vec),nrep)
for(l in 1:length(n_donor_vec)){
  mse_mat_overall = mse_mat_overall + mse_mat_list[[l]]
  #mse_mat_overall = mse_mat_overall + n_donor_vec[l] * mse_mat_list[[l]]
}
mse_mat_overall = mse_mat_overall/sum(n_donor_vec)
mse_mat_overall = sqrt(mse_mat_overall)


rbind(round(apply(mse_mat_overall, 1, mean),4),
      round(apply(mse_mat_overall, 1, sd),4)
)

mse_mat_list

pdf(paste('rmse_',nrep,'replications_bytask.pdf',sep = ""), width = 10, height = 5)
par(mgp = c(1,0.5,0), mfrow = c(2,4), mar = 0.5*c(2*0.9,0.2,2*1.95,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)
for(l in 1:length(n_donor_vec)){
  boxplot(t(sqrt(mse_mat_list[[l]]/n_donor_vec[l])), col = colvec, main = paste("Donor",l,'n=',n_donor_vec[l]), xaxt = "n", yaxt = "n", outline = F, cex.main = 1, lwd = 0.25)
  axis(2)
  legend('bottomright', legendvec, col = colvec, pch = 15, bty = "n", cex = 1.05)
}
dev.off()


mtlresearch_boxplot = boxplot(t(mse_mat_overall), plot = F)
ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))
par(mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)
boxplot(t(mse_mat_overall[1:length(fname_vec),]),col = colvec, main = "RMSE", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)
axis(2)
legend('bottomright', legendvec, col = colvec, pch = 15, bty = "n", cex = 1.05) 

pdf(paste('rmse_',nrep,'replications.pdf', sep = ""), width = 4, height = 3)
mtlresearch_boxplot = boxplot(t(mse_mat_overall), plot = F)
ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))

par(mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)

boxplot(t(mse_mat_overall[1:length(fname_vec),]),
        col = colvec, main = "RMSE", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)
axis(2)

legend('bottomright', legendvec,
       col = colvec, pch = 15, bty = "n", cex = 1.05) 
dev.off()









# 
# 
# mse_mat_overall_ = mse_mat_overall^2
# mtlresearch_boxplot = boxplot(t(mse_mat_overall_), plot = F)
# ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))
# par(mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)
# boxplot(t(mse_mat_overall_[1:length(fname_vec),]),col = colvec, main = "RMSE", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)
# axis(2)
# legend('bottomright', legendvec, col = colvec, pch = 15, bty = "n", cex = 1.05) 
# 
# pdf(paste('mse_',nrep,'replications.pdf', sep = ""), width = 4, height = 3)
# mse_mat_overall_ = mse_mat_overall^2
# mtlresearch_boxplot = boxplot(t(mse_mat_overall_), plot = F)
# ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))
# par(mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)
# boxplot(t(mse_mat_overall_[1:length(fname_vec),]),col = colvec, main = "RMSE", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)
# axis(2)
# legend('bottomright', legendvec, col = colvec, pch = 15, bty = "n", cex = 1.05) 
# dev.off()
