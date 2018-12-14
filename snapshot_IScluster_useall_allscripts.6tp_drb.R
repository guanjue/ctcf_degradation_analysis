scale_fc = function(x){
	xs = (log2(x))
	return((xs)-(xs[6]))
}

ave_fc = t(apply(ave+0.1, 1, scale_fc))

dr_nodrb = as.matrix(ave_fc[,-c(7:10)])

### snapshot
set.seed(2018)
not_used_col = -6
mod_all_bic <- densityMclust(as.vector(dr_nodrb[,not_used_col]), model='V')
set.seed(2018)
mod_all <- densityMclust(as.vector(dr_nodrb[,not_used_col]), G=3, model='V')

summary(mod_all)
attributes(mod_all)

cluster_id = mod_all$classification
cluster_mean = mod_all$parameters$mean
print('2nd GMM cluster means: ')
print(cluster_mean)

rainbow_cp = rev(rainbow(length(cluster_mean)))

rainbow_cp = c('blue', 'green', 'red')

pdf('gmm.density.pdf', width=14, height=7)
par(mfrow=c(1,2))
plot(mod_all, what = "density", data = as.vector(dr_nodrb[,not_used_col]), breaks = 50)
#plotDensityMclust1(mod_all, data = signal_mat_log2_vec_high, hist.col = "lightgrey", hist.border = "white",  breaks = "Sturges", type = "persp")

for (i in c(1:length(cluster_mean))){
	print(i)
	x_input = seq(-4,9, 0.1)
	cp_i_mean = mod_all$parameters$mean[i]
	cp_i_sd = (mod_all$parameters$variance$sigmasq[i])^0.5
	cp_i_pro = mod_all$parameters$pro[i]
	print(c(cp_i_mean, cp_i_sd, cp_i_pro))
	lines(x_input, cp_i_pro * dnorm(x_input, mean=cp_i_mean, sd=cp_i_sd), col=rainbow_cp[i])
}
plot(mod_all_bic, what = "BIC")
#plot(mod_all_bic, what = "density")
dev.off()

### get each cluster threshold
gmm_2nd_thresh = c()

f = function(x) dnorm(x, m=cluster_mean[1], sd=mod_all$parameters$variance$sigmasq[1]) * mod_all$parameters$pro[1] - dnorm(x, m=cluster_mean[2], sd=mod_all$parameters$variance$sigmasq[2]) * mod_all$parameters$pro[2]
gmm_2nd_thresh[1] = uniroot(f, interval=c(-1, 1))$root

f = function(x) dnorm(x, m=cluster_mean[3], sd=mod_all$parameters$variance$sigmasq[3]) * mod_all$parameters$pro[3] - dnorm(x, m=cluster_mean[2], sd=mod_all$parameters$variance$sigmasq[2]) * mod_all$parameters$pro[2]
gmm_2nd_thresh[2] = uniroot(f, interval=c(0, 5))$root


dr_nodrb_label = dr_nodrb[,not_used_col]
dr_nodrb_label[dr_nodrb_label>gmm_2nd_thresh[2]] = 10
dr_nodrb_label[(dr_nodrb_label<=gmm_2nd_thresh[2]) * (dr_nodrb_label>gmm_2nd_thresh[1]) >0] = 1
dr_nodrb_label[(dr_nodrb_label<=gmm_2nd_thresh[1])] = 0
dr_nodrb_label[dr_nodrb_label==10] = 2


### signal IS
ave_sig = ave[,1]
set.seed(2018)
mod_all_bic_sig <- densityMclust(log2(ave_sig+0.1), model='V')
set.seed(2018)
mod_all_sig <- densityMclust(log2(ave_sig+0.1), G=3, model='V')

cluster_id_sig = mod_all_sig$classification
cluster_mean_sig = mod_all_sig$parameters$mean
print('2nd GMM cluster means: ')
print(cluster_mean_sig)

rainbow_cp = rev(rainbow(length(cluster_mean_sig)))

rainbow_cp = c('blue', 'green', 'red')

pdf('gmm.density_sig.pdf', width=14, height=7)
par(mfrow=c(1,2))
plot(mod_all_sig, what = "density", data = log2(ave_sig+0.1), breaks = 50)
#plotDensityMclust1(mod_all, data = signal_mat_log2_vec_high, hist.col = "lightgrey", hist.border = "white",  breaks = "Sturges", type = "persp")

for (i in c(1:length(cluster_mean_sig))){
	print(i)
	x_input = seq(-4,9, 0.1)
	cp_i_mean = mod_all_sig$parameters$mean[i]
	cp_i_sd = (mod_all_sig$parameters$variance$sigmasq[i])^0.5
	cp_i_pro = mod_all_sig$parameters$pro[i]
	print(c(cp_i_mean, cp_i_sd, cp_i_pro))
	lines(x_input, cp_i_pro * dnorm(x_input, mean=cp_i_mean, sd=cp_i_sd), col=rainbow_cp[i])
}
plot(mod_all_bic_sig, what = "BIC")
#plot(mod_all_bic, what = "density")
dev.off()
gmm_2nd_thresh_sig = c()
f = function(x) dnorm(x, m=cluster_mean_sig[1], sd=mod_all_sig$parameters$variance$sigmasq[1]) * mod_all_sig$parameters$pro[1] - dnorm(x, m=cluster_mean_sig[2], sd=mod_all_sig$parameters$variance$sigmasq[2]) * mod_all_sig$parameters$pro[2]
gmm_2nd_thresh_sig[1] = uniroot(f, interval=c(0, 5))$root

f = function(x) dnorm(x, m=cluster_mean_sig[3], sd=mod_all_sig$parameters$variance$sigmasq[3]) * mod_all_sig$parameters$pro[3] - dnorm(x, m=cluster_mean_sig[2], sd=mod_all_sig$parameters$variance$sigmasq[2]) * mod_all_sig$parameters$pro[2]
gmm_2nd_thresh_sig[2] = uniroot(f, interval=c(0, 6))$root

ave_sig_label = log2(ave_sig+0.1)
ave_sig_label[log2(ave_sig+0.1)>gmm_2nd_thresh_sig[2]] = 10
ave_sig_label[(log2(ave_sig+0.1)<=gmm_2nd_thresh_sig[2]) * (log2(ave_sig+0.1)>gmm_2nd_thresh_sig[1]) >0] = 1
ave_sig_label[(log2(ave_sig+0.1)<=gmm_2nd_thresh_sig[1])] = 0
ave_sig_label[ave_sig_label==10] = 2


dr_nodrb_label_is = apply(cbind(ave_sig_label, dr_nodrb_label), 1, function(x) paste(x, collapse = '_'))
rare_thresh = 100

pdf('hist_IS.pdf')
hist(table(dr_nodrb_label_is), breaks=50)
abline(v=rare_thresh, col='red', lwd=1.5,lty=2)
dev.off()

#dr_nodrb_label_is_rare = rownames(table(dr_nodrb_label_is))[table(dr_nodrb_label_is)<rare_thresh]
dr_nodrb_label_is_rare = rownames(table(dr_nodrb_label_is)[table(dr_nodrb_label_is)<rare_thresh])
dr_nodrb_label_is_xxx = dr_nodrb_label_is
for (is_tmp in dr_nodrb_label_is_rare){
	dr_nodrb_label_is_xxx[dr_nodrb_label_is==is_tmp] = paste('x_x_x_x_x_', toString(i), sep='')
}
### qda
library(MASS)
qda_model = qda(cbind(log2(ave_sig+0.1), dr_nodrb[,not_used_col]), dr_nodrb_label_is_xxx)
new_IS = predict(qda_model,cbind(log2(ave_sig+0.1), dr_nodrb[,not_used_col]))$class
table_mat = cbind(table(dr_nodrb_label_is_xxx), table(new_IS))
#for (i in c(1:20)){
i=1
	print(i)
	qda_model = qda(cbind(log2(ave_sig+0.1), dr_nodrb[,not_used_col]), new_IS)
	new_IS = predict(qda_model,cbind(log2(ave_sig+0.1), dr_nodrb[,not_used_col]))$class
	table_mat = cbind(table_mat, table(new_IS))
#}

dr_nodrb_label_is_rare = rownames(table(new_IS)[table(new_IS)<rare_thresh])
dr_nodrb_label_is_xxx = new_IS

for (is_tmp in dr_nodrb_label_is_rare){
	dr_nodrb_label_is_xxx[new_IS==is_tmp] = paste('x_x_x_x_x_', toString(i), sep='')
}



plot_lim = 6
dr_IS_ave = log2(ave[order(dr_nodrb_label_is_xxx),]+1)

dr_IS_plot_ave = rbind(dr_IS_ave, ave_nopk[used_id_nopk,])

dr_IS_plot_ave[dr_IS_plot_ave>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('IS.s3.ave.', toString(nr), '.pdf', sep=''))
pheatmap(dr_IS_plot_ave, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()


dr_mclust_rep = log2(rep[order(dr_nodrb_label_is_xxx),]+1)
dr_IS_plot = rbind(dr_mclust_rep, rep_nopk[used_id_nopk,])

dr_IS_plot[dr_IS_plot>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('IS.s3.', toString(nr), '.pdf', sep=''))
pheatmap(dr_IS_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()


write.table(dr_nodrb_label_is_xxx, 'IS.classification.txt', quote=F, sep='\t', col.names=F, row.names=F)




signal_mat_index = signal_mat_log2
### get background index: '0'
signal_mat_index[signal_mat_log2<fdr_sig_thresh] = '0'
for (i in c(1:(length(cluster_mean)-1))){
	print(i)
	### get range id
	used_id_tmp = ( (signal_mat_log2>=gmm_2nd_thresh[i]) * (signal_mat_log2<gmm_2nd_thresh[i+1]) ) >0
	signal_mat_index[used_id_tmp] = toString(i)
}
### get top pk index
signal_mat_index[signal_mat_log2>=gmm_2nd_thresh[length(gmm_2nd_thresh)]] = toString(length(gmm_2nd_thresh))
### get index
signal_mat_index_vec = apply(signal_mat_index, 1, function(x) paste(x, collapse="_"))

write.table(signal_mat_index_vec, 'signal_mat_index_vec.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

### get index count
index_count = as.matrix(table(signal_mat_index_vec))
png('index_set_count_hist.png')
hist(log10(index_count))
dev.off()
index_count_enrich = as.matrix(index_count[index_count[,1]>100,])

write.table(index_count, 'index_count_all.txt', quote=FALSE, col.names=FALSE, row.names=TRUE, sep='\t')
write.table(index_count_enrich, 'index_count_enrich.txt', quote=FALSE, col.names=FALSE, row.names=TRUE, sep='\t')







