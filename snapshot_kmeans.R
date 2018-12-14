library(pheatmap)

rep1 = as.matrix(read.table('rep1_sig_TSnorm_mat.txt', header=TRUE))
rep2 = as.matrix(read.table('rep2_sig_TSnorm_mat.txt', header=TRUE))
ave = as.matrix(read.table('average_sig_TSnorm_mat.txt', header=TRUE))

rep1_nopk = as.matrix(read.table('rep1_sig_TSnorm_mat_nopk.txt', header=TRUE))
rep2_nopk = as.matrix(read.table('rep2_sig_TSnorm_mat_nopk.txt', header=TRUE))
ave_nopk = as.matrix(read.table('average_sig_TSnorm_mat_nopk.txt', header=TRUE))

rep1 = as.matrix(read.table('rep1_sig_s3norm_mat.txt', header=TRUE))
rep2 = as.matrix(read.table('rep2_sig_s3norm_mat.txt', header=TRUE))
ave = as.matrix(read.table('average_sig_s3norm_mat.txt', header=TRUE))

rep1_nopk = as.matrix(read.table('rep1_sig_s3norm_mat_nopk.txt', header=TRUE))
rep2_nopk = as.matrix(read.table('rep2_sig_s3norm_mat_nopk.txt', header=TRUE))
ave_nopk = as.matrix(read.table('average_sig_s3norm_mat_nopk.txt', header=TRUE))

set.seed(2018)
used_id_nopk = sample(dim(ave_nopk)[1], 5000)

pdf('signal_hist.pdf')
hist(ave, breaks=50)
dev.off()

pdf('signal_hist.log2.pdf')
hist(log2(ave), breaks=50)
dev.off()

tp_name = c('0hr', '4hr', '6hr', '12hr', '18hr', '24hr')

pdf('signal_hist.log2.all.pdf', height=14, width=7)
par(mfrow=c(6,1))
for (i in c(1:6)){
	sig = log2(ave[,i])
	breaks_vec = seq(min(sig)-1,max(sig)+1,length.out=10)
	hist(log2(ave[,i]), breaks=100, xlim=c(-5, 8), ylim=c(0,0.5), main=tp_name[i], freq=FALSE)
	box()
}
dev.off()


scale_fc = function(x){
	xs = (log2(x))
	return((xs)-(xs[6]))
}

ave_fc = t(apply(ave, 1, scale_fc))


pdf('fc_hist.pdf')
hist((ave_fc[,-1]), breaks=50)
dev.off()

pdf('fc_hist.log2.pdf')
hist(log2(ave_fc[,-dim(ave_fc)[2]]), breaks=100, xlim=c(-5.5, 15.5))
dev.off()

pdf('fc_hist.log2.all.pdf', height=14, width=7)
par(mfrow=c(5,1))
for (i in c(1:5)){
	hist((ave_fc[,i]), breaks=100, xlim=c(-4, 10), ylim=c(0,2.0), main=tp_name[i], freq=FALSE)
	abline(v=0, col='red',lwd=1.5)
	box()
}
dev.off()


kmeans_dist = c()
ave_fc_col_var = sum(apply(ave_fc, 2, var))

for (i in c(1:20)){
print(i)
#nr = 10
set.seed(2018)
dr = as.matrix(ave_fc)
fit = kmeans(dr, i)

dr_tmp_dist_all = c()
for (j in c(1:i)){
	print(j)
	dr_tmp = dr[fit$cluster==j,]
	dr_tmp_mean = colMeans(dr_tmp)
	dr_tmp_dist = t(apply(dr_tmp, 1, function(x) (x-dr_tmp_mean)^2))
	dr_tmp_dist_all = rbind(dr_tmp_dist_all, dr_tmp_dist)
}
print(dim(dr_tmp_dist_all))
dr_tmp_dist_all_sum = sum(dr_tmp_dist_all)/dim(dr_tmp_dist_all)[1]/ave_fc_col_var
kmeans_dist[i] = dr_tmp_dist_all_sum
}

pdf('kmeans_dist.pdf')
plot(c(1:length(kmeans_dist)), kmeans_dist, lwd = 1.5, type='l')
points(c(1:length(kmeans_dist)), kmeans_dist, col='black')
dev.off()





set.seed(2018)
nr = 15
dr = as.matrix(ave_fc)
fit = kmeans(dr, nr)
sig_reps = c()
sig_reps_nopk = c()
for (i in c(1:dim(rep1)[2])){
	sig_reps = cbind(sig_reps, rep1[,i], rep2[,i])
	sig_reps_nopk = cbind(sig_reps_nopk, rep1_nopk[,i], rep2_nopk[,i])
}


dr_kmeans = log2(sig_reps[order(fit$cluster),]+1)

my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.s3.', toString(nr), '.pdf', sep=''))
pheatmap(rbind(dr_kmeans, sig_reps_nopk[used_id_nopk,]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

dr_kmeans_ave = log2(ave[order(fit$cluster),]+1)

my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.s3.ave.', toString(nr), '.pdf', sep=''))
pheatmap(rbind(dr_kmeans_ave, ave_nopk[used_id_nopk,]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

png(paste('kmean.s3.ave.', toString(nr), '.png', sep=''))
pheatmap(rbind(dr_kmeans_ave, ave_nopk[used_id_nopk,]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()






my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.TS.', toString(nr), '.pdf', sep=''))
pheatmap(rbind(dr_kmeans, sig_reps_nopk[used_id_nopk,]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

dr_kmeans_ave = log2(ave[order(fit$cluster),]+1)

my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.TS.ave.', toString(nr), '.pdf', sep=''))
pheatmap(rbind(dr_kmeans_ave, ave_nopk[used_id_nopk,]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

png(paste('kmean.TS.ave.', toString(nr), '.png', sep=''))
pheatmap(rbind(dr_kmeans_ave, ave_nopk[used_id_nopk,]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()


