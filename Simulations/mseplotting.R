
methodvec = c('MTL_transfer','MTL_transfer_debiased',
              'MTL_notransfer','MTL_notransfer_merged',
              'STL_transfer','STL_transfer_debiased') #,'STL_notransfer','STL_notransfer_merged'

msemat = matrix(NA, nrow = length(methodvec), ncol = 10) ; p = 100 ; K = 50

for(i in 1:length(methodvec)){
  load(paste(methodvec[i], '/summary.Rdata',sep=""))
  msemat[i,] = log(p*K*estimation_msevec^2)
}

colvec = c('red1','pink2', 'green1','green4', 'blue','skyblue3')

pdf('logfrobnorm.pdf', width = 5.75, height = 3.5)
mtlresearch_boxplot = boxplot(t(msemat), plot = F)
ylim_mtlresearch = c(min(mtlresearch_boxplot$stats), max(mtlresearch_boxplot$stats))

par(mgp = c(1,0.5,0), mar = 0.5*c(2*0.65,0.2,2*1.75,0.2), oma = 0.5*c(0.1, 5.75, 1, 0.1), cex.main = 1.15, cex.axis = 1.15)

boxplot(t(msemat),
        col = colvec, main = "log Frobenius-norm", xaxt = "n", yaxt = "n", outline = F, cex.main = 1.35, ylim = ylim_mtlresearch)
abline(v = 2.5, lty = "longdash") ; abline(v = 4.5, lty = "longdash") ; #abline(v = 6.5, lty = "longdash")
axis(2)
legend('topleft', c('TMTL(Fused)','TMTL(Debiased)',
                    'MTL(Target)','MTL(Full)',
                    'TSTL(Fused)','TSTL(Debiased)'),
       col = colvec, pch = 15, bty = "n", cex = 1) 

dev.off()