###########
R2 = function(d_IP_ts, qPCR_expand){
        r2 = 1-sum((qPCR_expand-d_IP_ts)^2) / sum((qPCR_expand-mean(qPCR_expand))^2)
        return(r2)
}
###########

a = c(1:10)
b = c(1:10)
b2 = b*2-1

pdf('R2_vs_R.pdf', width=6, height=3.5)
par(mfrow=c(1,2))
plot(a,b, xlim=c(1,20), ylim=c(0,20), main=paste('R2: ', round(R2(a,b),2), '; ', 'Correlation: ', cor(a,b), sep=''))
abline(0,1,col='red')

plot(a,b2, xlim=c(1,20), ylim=c(0,20), main=paste('R2: ', round(R2(a,b2),2), '; ', 'Correlation: ', cor(a,b), sep=''))
abline(0,1,col='red')
dev.off()