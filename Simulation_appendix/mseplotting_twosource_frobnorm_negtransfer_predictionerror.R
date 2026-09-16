
settingvec = c('twosource_balanced','twosource_aligned')    

methodvec = c('MTL_transfer','MTL_transfer_debiased',
              'MTL_notransfer','MTL_notransfer_merged',
              'STL_transfer','STL_transfer_debiased')              

msemat = matrix(NA, nrow = length(methodvec)*length(settingvec), ncol = 20) ; p = 100 ; K = 50
predmat = matrix(NA, nrow = length(methodvec)*length(settingvec), ncol = 20)          

for(i in 1:length(settingvec)){
  for(j in 1:length(methodvec)){
    #cat('i = ', i, ', j = ',j,'\n')
    
    load(paste(settingvec[i],'/',methodvec[j],'/summary.Rdata', sep = ""))
    msemat[length(methodvec)*(i-1)+j, 1:20] = log(p*estimation_msevec^2)               
    predmat[length(methodvec)*(i-1)+j, 1:20] = prediction_msevec                       
    
  }
}

colvec = c('red1','pink2', 'green1','green4', 'blue','skyblue3')




#### negative transfer frequency (vs MTL(Target), paired by replication) ####         
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
cat('\n--- negative transfer frequency (prediction RMSE) ---\n') ; print(round(ntmat_pred, 2))




## Balanced 
apply(t(msemat[0*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[0*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)

## aligned
apply(t(msemat[1*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[1*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)


#### prediction error summary ####                                                    
## Balanced 
apply(t(predmat[0*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(4)
apply(t(predmat[0*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(4)

## aligned
apply(t(predmat[1*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(4)
apply(t(predmat[1*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(4)






pdf('logfrobnorm_twosource.pdf', width = length(settingvec)*2.5, height = 3)
mtlresearch_boxplot = boxplot(t(msemat), plot = F)
ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))

par(mfrow = c(1,length(settingvec)), mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)

boxplot(t(msemat[1:length(methodvec),]),
        col = colvec, main = "Balanced shifts", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)
axis(2)

boxplot(t(msemat[1*length(methodvec)+1:length(methodvec),]),
        col = colvec, main = "Aligned shifts", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)

legend('bottomright', c('TMTL(Fused)','TMTL(Debiased)',
                        'MTL(Target)','MTL(Full)',
                        'TSTL(Fused)','TSTL(Debiased)'),
       col = colvec, pch = 15, bty = "n", cex = 0.85) 

dev.off()


#### prediction error plot (same layout) ####                                         
pdf('predrmse_twosource.pdf', width = length(settingvec)*2.5, height = 3)
pred_boxplot = boxplot(t(predmat), plot = F)
ylim_pred = c(min(pred_boxplot$stats), max(pred_boxplot$stats))

par(mfrow = c(1,length(settingvec)), mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)

boxplot(t(predmat[1:length(methodvec),]),
        col = colvec, main = "Balanced shifts", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_pred, lwd = 0.25)
axis(2)

boxplot(t(predmat[1*length(methodvec)+1:length(methodvec),]),
        col = colvec, main = "Aligned shifts", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_pred, lwd = 0.25)

legend('bottomright', c('TMTL(Fused)','TMTL(Debiased)',
                        'MTL(Target)','MTL(Full)',
                        'TSTL(Fused)','TSTL(Debiased)'),
       col = colvec, pch = 15, bty = "n", cex = 0.85) 

dev.off()