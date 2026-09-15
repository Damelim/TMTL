library(dplyr)
settingvec = c('threesource_noskew','threesource_aligned_shift0.25','threesource_aligned_shift0.5','threesource_aligned_shift1')              #,'twosource_verybigskew')

methodvec = c('MTL_transfer','MTL_transfer_debiased',
              'MTL_notransfer','MTL_notransfer_merged',
              'STL_transfer_cluster_sparse','STL_transfer_debiased')             #,'STL_notransfer','STL_notransfer_merged'

msemat = matrix(NA, nrow = length(methodvec)*length(settingvec), ncol = 100) ; p = 200 ; K = 20
predmat = matrix(NA, nrow = length(methodvec)*length(settingvec), ncol = 100)          # <-- added

#msemat = data.frame(msemat)

for(i in 1:length(settingvec)){
  for(j in 1:length(methodvec)){
    #rownames(msemat)[length(methodvec)*(i-1)+j] = paste('setting',settingvec[i],'method',methodvec[j])
    load(paste(settingvec[i],'/',methodvec[j],'/summary.Rdata', sep = ""))
    msemat[length(methodvec)*(i-1)+j, 1:100] = log(p*estimation_msevec[1:100]^2)          # <-- p*K -> p
    predmat[length(methodvec)*(i-1)+j, 1:100] = prediction_msevec[1:100]                  # <-- added
  }
}


print(sum(is.na(msemat)))




#### negative transfer frequency (vs MTL(Target), paired by replication) ####         # <-- added
ntmat_est = matrix(NA, nrow = length(settingvec), ncol = length(methodvec))
ntmat_pred = matrix(NA, nrow = length(settingvec), ncol = length(methodvec))
rownames(ntmat_est) = settingvec ; colnames(ntmat_est) = methodvec
rownames(ntmat_pred) = settingvec ; colnames(ntmat_pred) = methodvec
for(i in 1:length(settingvec)){
  refrow = length(methodvec)*(i-1) + which(methodvec == 'MTL_notransfer')
  for(j in 1:length(methodvec)){
    ntmat_est[i,j] = mean(msemat[length(methodvec)*(i-1)+j, ] > msemat[refrow, ])
    ntmat_pred[i,j] = mean(predmat[length(methodvec)*(i-1)+j, ] > predmat[refrow, ])
  }
}
cat('\n--- negative transfer frequency (coefficient error) ---\n') ; print(round(ntmat_est, 2))
#cat('\n--- negative transfer frequency (prediction RMSE) ---\n') ; print(round(ntmat_pred, 2))




## Balanced 
apply(t(msemat[0*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[0*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)

## alpha0.25
apply(t(msemat[1*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[1*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)

## alpha0.5
apply(t(msemat[2*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[2*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)

## alpha1
apply(t(msemat[3*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[3*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)


#### prediction error summary ####                                                    # <-- added
## Balanced 
apply(t(predmat[0*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(predmat[0*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)

## alpha0.25
apply(t(predmat[1*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(predmat[1*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)

## alpha0.5
apply(t(predmat[2*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(predmat[2*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)

## alpha1
apply(t(predmat[3*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(predmat[3*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)






colvec = c('red1','pink2', 'green1','green4', 'blue','skyblue3')

pdf('logfrobnorm_threesource.pdf', width = length(settingvec)*2.5, height = 3)
mtlresearch_boxplot = boxplot(t(msemat), plot = F)
ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))

par(mfrow = c(1,length(settingvec)), mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)

boxplot(t(msemat[1:length(methodvec),]),
        col = colvec, main = "Balanced shifts", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)
axis(2)

boxplot(t(msemat[1*length(methodvec)+1:length(methodvec),]),
        col = colvec, main = expression(paste(alpha,'=0.25')), xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)

boxplot(t(msemat[2*length(methodvec)+1:length(methodvec),]),
        col = colvec, main = expression(paste(alpha,'=0.5')), xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)

boxplot(t(msemat[3*length(methodvec)+1:length(methodvec),]),
        col = colvec, main = expression(paste(alpha,'=1')), xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)

legend('bottomright', c('TMTL(Fused)','TMTL(Debiased)',
                        'MTL(Target)','MTL(Full)',
                        'TSTL(Fused)','TSTL(Debiased)'),
       col = colvec, pch = 15, bty = "n", cex = 0.85) 

dev.off()


#### prediction error plot (same layout) ####                                         # <-- added
pdf('predrmse_threesource.pdf', width = length(settingvec)*2.5, height = 3)
pred_boxplot = boxplot(t(predmat), plot = F)
ylim_pred = c(min(pred_boxplot$stats), max(pred_boxplot$stats))

par(mfrow = c(1,length(settingvec)), mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)

boxplot(t(predmat[1:length(methodvec),]),
        col = colvec, main = "Balanced shifts", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_pred, lwd = 0.25)
axis(2)

boxplot(t(predmat[1*length(methodvec)+1:length(methodvec),]),
        col = colvec, main = expression(paste(alpha,'=0.25')), xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_pred, lwd = 0.25)

boxplot(t(predmat[2*length(methodvec)+1:length(methodvec),]),
        col = colvec, main = expression(paste(alpha,'=0.5')), xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_pred, lwd = 0.25)

boxplot(t(predmat[3*length(methodvec)+1:length(methodvec),]),
        col = colvec, main = expression(paste(alpha,'=1')), xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_pred, lwd = 0.25)

legend('bottomright', c('TMTL(Fused)','TMTL(Debiased)',
                        'MTL(Target)','MTL(Full)',
                        'TSTL(Fused)','TSTL(Debiased)'),
       col = colvec, pch = 15, bty = "n", cex = 0.85) 

dev.off()