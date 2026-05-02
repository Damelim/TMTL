# settingvec = c('alpha0','alpha1:5','alpha1:2','alpha1')
# 
# methodvec = c('MTL_transfer','MTL_transfer_debiased',
#               'MTL_notransfer','MTL_notransfer_merged',
#               'STL_transfer','STL_transfer_debiased') #,'STL_notransfer','STL_notransfer_merged'
# 
# msemat = matrix(NA, nrow = length(methodvec)*length(settingvec), ncol = 10) ; p = 100 ; K = 50
# 
# 
# for(i in 1:length(settingvec)){
#   for(j in 1:length(methodvec)){
#     cat('i = ', i, ', j = ',j,'\n')
#     load(paste(settingvec[i],'/',methodvec[j],'/summary.Rdata', sep = ""))
#     msemat[length(methodvec)*(i-1)+j, ] = log(p*K*estimation_msevec^2)
#   }
# }
# 
# colvec = c('red1','pink2', 'green1','green4', 'blue','skyblue3')
# 
# pdf('logfrobnorm_onesource.pdf', width = length(settingvec)*2, height = 3)
# mtlresearch_boxplot = boxplot(t(msemat), plot = F)
# ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))
# 
# par(mfrow = c(1,length(settingvec)), mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 1, 0.1), cex.main = 1.5, cex.axis = 1.5)
# 
# boxplot(t(msemat[1:length(methodvec),]),
#         col = colvec, main = expression(paste(alpha,"=0")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_mtlresearch, lwd = 0.25)
# #abline(v = 2.5, lty = "longdash") ; abline(v = 4.5, lty = "longdash") ; #abline(v = 6.5, lty = "longdash")
# axis(2)
# legend('topleft', c('TMTL(Fused)','TMTL(Debiased)',
#                     'MTL(Target)','MTL(Full)',
#                     'TSTL(Fused)','TSTL(Debiased)'),
#        col = colvec, pch = 15, bty = "n", cex = 1.25) 
# 
# boxplot(t(msemat[1*length(methodvec)+1:length(methodvec),]),
#         col = colvec, main = expression(paste(alpha,"=1/5")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_mtlresearch, lwd = 0.25)
# 
# boxplot(t(msemat[2*length(methodvec)+1:length(methodvec),]),
#         col = colvec, main = expression(paste(alpha,"=1/2")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_mtlresearch, lwd = 0.25)
# 
# boxplot(t(msemat[3*length(methodvec)+1:length(methodvec),]),
#         col = colvec, main = expression(paste(alpha,"=1")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_mtlresearch, lwd = 0.25)
# 
# dev.off()

















settingvec = c('alpha0','alpha1:5','alpha1:2','alpha1')

methodvec = c('MTL_transfer','MTL_transfer_debiased',
              'MTL_notransfer','MTL_notransfer_merged',
              'STL_transfer','STL_transfer_debiased') #,'STL_notransfer','STL_notransfer_merged'

msemat = matrix(NA, nrow = length(methodvec)*length(settingvec), ncol = 20) ; p = 100 ; K = 50


for(i in 1:length(settingvec)){
  for(j in 1:length(methodvec)){
    cat('i = ', i, ', j = ',j,'\n')
    load(paste(settingvec[i],'/',methodvec[j],'/summary.Rdata', sep = ""))
    msemat[length(methodvec)*(i-1)+j, 1:10] = log(p*K*estimation_msevec^2)
    
    load(paste(settingvec[i],'_diffseed/',methodvec[j],'/summary.Rdata', sep = ""))
    msemat[length(methodvec)*(i-1)+j, 11:20] = log(p*K*estimation_msevec[11:20]^2)
  }
}

colvec = c('red1','pink2', 'green1','green4', 'blue','skyblue3')

pdf('logfrobnorm_onesource.pdf', width = length(settingvec)*2, height = 3)
mtlresearch_boxplot = boxplot(t(msemat), plot = F)
ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))

par(mfrow = c(1,length(settingvec)), mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.25, 3.75, 1, 0.1), cex.main = 1.5, cex.axis = 1.5)

boxplot(t(msemat[1:length(methodvec),]),
        col = colvec, main = expression(paste(alpha,"=0")), xaxt = "n", yaxt = "n", outline = F, cex.main = 1.5, ylim = ylim_mtlresearch, lwd = 0.25)
#abline(v = 2.5, lty = "longdash") ; abline(v = 4.5, lty = "longdash") ; #abline(v = 6.5, lty = "longdash")
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