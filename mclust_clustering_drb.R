library(pheatmap)



rep0 = (read.table('DPGP_TSnorm/rep_sig_TSnorm_mat.txt.bed', header=TRUE))

rep = as.matrix(rep0[,-c(1:4, 7:18)])



rep_mean_old = function(x){
	xm = c((x[1]+x[2])/2, (x[3]+x[4])/2, (x[5]+x[6])/2, (x[7]+x[8])/2, (x[11]+x[10]+x[11])/3)	
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
	return((xs)/(xs[1]))
}


scale_fc_log2 = function(x){
	xs = (log2(x))
	return(((xs)-(xs[1]))[-1])
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

rep_mean_fc = t(apply(rep_mean+small_num,1,scale_fc_log2))


library(mclust)

set.seed(2018)
nr = 10

plot_lim = 8

BIC = mclustBIC(rep_mean_fc, modelNames=c("EII", "VII"))
png('bic.drb.png')
plot(BIC)
dev.off()

pdf('bic.drb.pdf')
plot(BIC)
dev.off()


fit = Mclust(rep_mean_fc, x = BIC)

nr1= dim(table(fit$classification))

fit_mean = c()
for (i in c(1:nr1)){
	sig_tmp = rep[fit$classification==i,]
	sig_tmp_mean = mean(sig_tmp)
	fit_mean[i] = sig_tmp_mean
}



plot_cluster_rank = order(fit_mean)

dr_mclust_plot = c()
dr_mclust_plot_label = c()
dr_mclust_plot_cluster_num = c()
pk_k_all = c()

k=0
for (i in plot_cluster_rank){
	k = k+1
	dr_mclust_plot_cluster_num = c(dr_mclust_plot_cluster_num, rep(k, sum(fit$classification==i)))
	print(sum(fit$classification==i))
	dr_mclust_plot = rbind(dr_mclust_plot, log2(rep[fit$classification==i,]+1))
	dr_mclust_plot_label = rbind(dr_mclust_plot_label, cbind(rep(i, sum(fit$classification==i)),rep(i, sum(fit$classification==i))) )

	pk_k = rep0[fit$classification==i,]
	pk_k_all = rbind(pk_k_all, cbind(pk_k, rep(k, sum(fit$classification==i))))
	write.table(pk_k, paste('mclust.TS.drb.', toString(k), '.pk.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
	### plot boxplot
	pdf(paste('mclust.box.drb.', toString(k), '.pdf', sep=''))
	sigmat0 = as.matrix((rep[fit$classification==i,]+small_num))
	sigmat = t(apply(sigmat0, 1, scale_fc_log2_mean))
	boxplot(sigmat, ylim=c(-3, 1) ,outline=FALSE)
	if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	lines(c(1:(dim(sigmat)[2])), colMedian(sigmat), col='black')
	dev.off()
}



dr_mclust_plot[dr_mclust_plot>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.TS.drb.nodrb.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

pdf(paste('mclust.TS.drb.nodrb.label.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot_label, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

png(paste('mclust.TS.drb.nodrb.label.', toString(nr), '.png', sep=''))
pheatmap(dr_mclust_plot_label, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

colnames(pk_k_all)[length(colnames(pk_k_all))] = 'cluster_id'
write.table(pk_k_all, paste('mclust.TS.drb.all.pk.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')



