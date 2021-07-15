#setwd('/path to your working directory/')

sum5s_avg =read.table('sum5s_avg')
sum5s_avg=as.matrix(sum5s_avg, ncol=2, byrow=580)
sum3s_avg=read.table('sum3s_avg')
sum3s_avg=as.matrix(sum3s_avg, ncol=2, byrow=580)

plot(sum5s_avg[,1], sum5s_avg[,2], type='l', cex=0.5, col='red', lwd=1, ylab='Nucleosome occupancy (mean count)', lty=1, xlab="Relative positions with respect to internal exons",xaxt='n', xlim=c(-610,610), ylim=c(0,1.2))
points (sum3s_avg[,1], sum3s_avg[,2], type='l', pch=20, cex=0.5, col='red', lwd=1, lty=1, ylab='count_avg', xlab='bp')

legend(320, 1.2, legend=c("H3K4me3"), col=c("red"), lty=c(1), cex=0.8)
legend(-580,1.2, legend=bquote(LncRNA_exons), 
       bty="n", inset=c(0, .01), text.font=1.8, cex = 1.2)

ticks = c(-580,-80, 130, 610)
axis(1, at=ticks, labels=c('580','0', '0', '580'))

abline(v = 130, col="black", lwd=0.5, lty=2)
abline(v = -80, col="black", lwd=0.5, lty=2)

