library(pheatmap)



rep0 = (read.table('DPGP_TSnorm/rep_sig_TSnorm_mat.txt.bed', header=TRUE))

rep = as.matrix(rep0[,-c(1:4, 7:18)])



rep_mean_old = function(x){
	xm = c((x[1]+x[2])/2, (x[3]+x[4])/2, (x[5]+x[6])/2, (x[7]+x[8])/2, (x[11]+x[10])/2)	
	return(xm)
}


rep_mean_fun = function(x){
	xm = c()
	for (i in c(1:(length(x)/2))){
		xm_tmp = (x[1+(i-1)*2]+x[2+(i-1)*2])/2
		xm[i] = xm_tmp
	}	
	return(xm)
}


scale_fc_mean = function(x){
	xs = ((x))
	return((xs)/(xs[1]/2+xs[2]/2))
}

scale_fc = function(x){
	xs = ((x))
	return(((xs)/(xs[1]))[-1])
}


scale_fc_log2 = function(x){
	xs = (log2(x))
	return((xs)-(xs[1]))
}


scale_fc_log2_mean = function(x){
	xs = ((x))
	return(log2((xs)/(xs[1]/2+xs[2]/2)))
}


colMedian = function(x){
	xm = apply(x,2,median)
	return(xm)
}


small_num = 10
rep_mean = t(apply(rep, 1, rep_mean_old))

rep_mean_fc = t(apply(rep_mean+small_num,1,scale_fc))

kmeans_dist = c()
ave_fc_col_var = sum(apply(rep_mean_fc, 2, var))

for (i in c(1:20)){
print(i)
#nr = 10
set.seed(2018)
dr = as.matrix(rep_mean_fc)
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

pdf('kmeans_dist.drb.pdf')
plot(c(1:length(kmeans_dist)), kmeans_dist, lwd = 1.5, type='l')
points(c(1:length(kmeans_dist)), kmeans_dist, col='black')
dev.off()


set.seed(2018)
nr = 10
nr1= 12

plot_lim = 8
### cluster without drb
fit = kmeans(rep_mean_fc, nr1)

fit_mean = c()
for (i in c(1:nr1)){
	sig_tmp = rep[fit$cluster==i,]
	sig_tmp_mean = mean(sig_tmp[,])
	fit_mean[i] = sig_tmp_mean
}


#pdf('hist_rep.pdf')
#hist(log10(rep), breaks=50)
#dev.off()

plot_cluster_rank = order(fit_mean)

dr_kmeans_plot = c()
dr_kmeans_plot_label = c()
dr_kmeans_plot_cluster_num = c()
pk_k_all = c()

k=0
for (i in plot_cluster_rank){
	k = k+1
	dr_kmeans_plot_cluster_num = c(dr_kmeans_plot_cluster_num, rep(k, sum(fit$cluster==i)))
	print(sum(fit$cluster==i))
	dr_kmeans_plot = rbind(dr_kmeans_plot, log2(rep[fit$cluster==i,]+1))
	dr_kmeans_plot_label = rbind(dr_kmeans_plot_label, cbind(rep(i, sum(fit$cluster==i)),rep(i, sum(fit$cluster==i))) )

	pk_k = rep0[fit$cluster==i,]
	pk_k_all = rbind(pk_k_all, cbind(pk_k, rep(k, sum(fit$cluster==i))))
	write.table(pk_k, paste('kmean.TS.drb.', toString(k), '.pk.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
	### plot boxplot
	pdf(paste('kmeans.drb.box.', toString(k), '.pdf', sep=''))
	sigmat0 = as.matrix((rep[fit$cluster==i,]+small_num))
	sigmat = t(apply(sigmat0, 1, scale_fc_log2_mean))
	boxplot(sigmat, ylim=c(-3, 1) ,outline=FALSE)
	if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	lines(c(1:(dim(sigmat)[2])), colMedian(sigmat), col='black')
	dev.off()
}

dr_kmeans_plot[dr_kmeans_plot>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.TS.drb.', toString(nr), '.pdf', sep=''))
pheatmap(dr_kmeans_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

pdf(paste('kmean.TS.drb.label.', toString(nr), '.pdf', sep=''))
pheatmap(dr_kmeans_plot_label, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

png(paste('kmean.TS.drb.label.', toString(nr), '.png', sep=''))
pheatmap(dr_kmeans_plot_label, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

colnames(pk_k_all)[length(colnames(pk_k_all))] = 'cluster_id'
write.table(pk_k_all, paste('kmean.TS.drb.all.pk.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')










