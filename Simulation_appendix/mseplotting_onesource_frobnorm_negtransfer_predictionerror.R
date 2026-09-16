settingvec = c('alpha0','alpha1:5','alpha1:2','alpha1')
methodvec = c('MTL_transfer','MTL_transfer_debiased',
              'MTL_notransfer','MTL_notransfer_merged',
              'STL_transfer','STL_transfer_debiased')  
msemat = matrix(NA, nrow = length(methodvec)*length(settingvec), ncol = 20) ; p = 100 ; K = 50
predmat = matrix(NA, nrow = length(methodvec)*length(settingvec), ncol = 20)   
for(i in 1:length(settingvec)){
  for(j in 1:length(methodvec)){

    if(j %in% 1:4){
      
      load(paste(settingvec[i],'/',methodvec[j],'/summary.Rdata', sep = ""))
      msemat[length(methodvec)*(i-1)+j, 1:20] = log(p*estimation_msevec^2)         
      predmat[length(methodvec)*(i-1)+j, 1:20] = prediction_msevec                 
      
    }else{
      
      load(paste(settingvec[i],'/',methodvec[j],'/summary.Rdata', sep = ""))
      msemat[length(methodvec)*(i-1)+j, 1:20] = log(p*estimation_msevec^2)          
      predmat[length(methodvec)*(i-1)+j, 1:20] = prediction_msevec                 
    }
    
  } # end j 
} # end i


library(dplyr)

colvec = c('red1','pink2', 'green1','green4', 'blue','skyblue3')
if(sum(is.na(msemat)) == 0){
  pdf('logfrobnorm_onesource.pdf', width = length(settingvec)*2, height = 3)
  mtlresearch_boxplot = boxplot(t(msemat), plot = F)
  ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))
  
  par(mfrow = c(1,length(settingvec)), mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 1, 0.1), cex.main = 1.5, cex.axis = 1.5)
  
  boxplot(t(msemat[1:length(methodvec),]),
          col = colvec, main = expression(paste(alpha,"=0")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_mtlresearch, lwd = 0.25)
  axis(2)
  legend('topleft', c('TMTL(Fused)','TMTL(Debiased)',
                      'MTL(Target)','MTL(Full)',
                      'TSTL(Fused)','TSTL(Debiased)'),
         col = colvec, pch = 15, bty = "n", cex = 1.25) 
  
  boxplot(t(msemat[1*length(methodvec)+1:length(methodvec),]),
          col = colvec, main = expression(paste(alpha,"=1/5")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_mtlresearch, lwd = 0.25)
  
  boxplot(t(msemat[2*length(methodvec)+1:length(methodvec),]),
          col = colvec, main = expression(paste(alpha,"=1/2")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_mtlresearch, lwd = 0.25)
  
  boxplot(t(msemat[3*length(methodvec)+1:length(methodvec),]),
          col = colvec, main = expression(paste(alpha,"=1")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_mtlresearch, lwd = 0.25)
  
  dev.off()
}else{
  print('Error')
}

#### prediction error plot (same layout) ####                                    
if(sum(is.na(predmat)) == 0){
  pdf('predrmse_onesource.pdf', width = length(settingvec)*2, height = 3)
  pred_boxplot = boxplot(t(predmat), plot = F)
  ylim_pred = c(min(pred_boxplot$stats), max(pred_boxplot$stats))
  
  par(mfrow = c(1,length(settingvec)), mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 1, 0.1), cex.main = 1.5, cex.axis = 1.5)
  
  boxplot(t(predmat[1:length(methodvec),]),
          col = colvec, main = expression(paste(alpha,"=0")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_pred, lwd = 0.25)
  axis(2)
  legend('topleft', c('TMTL(Fused)','TMTL(Debiased)',
                      'MTL(Target)','MTL(Full)',
                      'TSTL(Fused)','TSTL(Debiased)'),
         col = colvec, pch = 15, bty = "n", cex = 1.25) 
  
  boxplot(t(predmat[1*length(methodvec)+1:length(methodvec),]),
          col = colvec, main = expression(paste(alpha,"=1/5")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_pred, lwd = 0.25)
  
  boxplot(t(predmat[2*length(methodvec)+1:length(methodvec),]),
          col = colvec, main = expression(paste(alpha,"=1/2")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_pred, lwd = 0.25)
  
  boxplot(t(predmat[3*length(methodvec)+1:length(methodvec),]),
          col = colvec, main = expression(paste(alpha,"=1")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_pred, lwd = 0.25)
  
  dev.off()
}else{
  print('Error')
}




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






## alpha = 0
apply(t(msemat[0*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[0*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)
## alpha = 1/5
apply(t(msemat[1*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[1*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)
## alpha = 1/2
apply(t(msemat[2*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[2*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)
## alpha = 11
apply(t(msemat[3*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(3)
apply(t(msemat[3*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(3)




#### prediction error summary ####                                              
## alpha = 0
apply(t(predmat[0*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(4)
apply(t(predmat[0*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(4)
## alpha = 1/5
apply(t(predmat[1*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(4)
apply(t(predmat[1*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(4)
## alpha = 1/2
apply(t(predmat[2*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(4)
apply(t(predmat[2*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(4)
## alpha = 1
apply(t(predmat[3*length(methodvec)+1:length(methodvec),]),2,mean) %>% round(4)
apply(t(predmat[3*length(methodvec)+1:length(methodvec),]),2,sd) %>% round(4)