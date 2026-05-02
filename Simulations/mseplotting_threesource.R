settingvec = c('threesource_noskew','threesource_bigskew')              #,'twosource_verybigskew')

methodvec = c('MTL_transfer','MTL_transfer_debiased',
              'MTL_notransfer','MTL_notransfer_merged',
              'STL_transfer','STL_transfer_debiased')             #,'STL_notransfer','STL_notransfer_merged'

msemat = matrix(NA, nrow = length(methodvec)*length(settingvec), ncol = 20) ; p = 100 ; K = 50


for(i in 1:length(settingvec)){
  for(j in 1:length(methodvec)){
    cat('i = ', i, ', j = ',j,'\n')
    load(paste(settingvec[i],'/',methodvec[j],'/summary.Rdata', sep = ""))
    msemat[length(methodvec)*(i-1)+j, 1:3] = log(p*K*estimation_msevec[1:3]^2)
    
    load(paste(settingvec[i],'_diffseed/',methodvec[j],'/summary.Rdata', sep = ""))
    msemat[length(methodvec)*(i-1)+j, 4:10] = log(p*K*estimation_msevec[4:10]^2)
    
    load(paste(settingvec[i],'_diffseed2/',methodvec[j],'/summary.Rdata', sep = ""))
    msemat[length(methodvec)*(i-1)+j, 11:15] = log(p*K*estimation_msevec[11:15]^2)
    
    load(paste(settingvec[i],'_diffseed3/',methodvec[j],'/summary.Rdata', sep = ""))
    msemat[length(methodvec)*(i-1)+j, 16:20] = log(p*K*estimation_msevec[16:20]^2)
  }
}

colvec = c('red1','pink2', 'green1','green4', 'blue','skyblue3')

pdf('logfrobnorm_threesource.pdf', width = length(settingvec)*2.5, height = 3)
mtlresearch_boxplot = boxplot(t(msemat), plot = F)
ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))

par(mfrow = c(1,length(settingvec)), mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 0.25, 0.1), cex.main = 1.05, cex.axis = 1.05)

boxplot(t(msemat[1:length(methodvec),]),
        col = colvec, main = "No skew", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)
axis(2)

boxplot(t(msemat[1*length(methodvec)+1:length(methodvec),]),
        col = colvec, main = "Big skew", xaxt = "n", yaxt = "n", outline = F, cex.main = 1, ylim = ylim_mtlresearch, lwd = 0.25)

legend('bottomright', c('TMTL(Fused)','TMTL(Debiased)',
                        'MTL(Target)','MTL(Full)',
                        'TSTL(Fused)','TSTL(Debiased)'),
       col = colvec, pch = 15, bty = "n", cex = 0.85) 

dev.off()
