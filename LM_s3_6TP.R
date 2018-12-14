library(mclust)
library(pheatmap)

#######################################
### FOR 6 time points
#######################################
print('FOR 6 time points')
all_sample_input_table = read.table('bam_bw_bedgraph/file_list.all.txt', header=F)[c(1:15),]
IP_list = all_sample_input_table[,1]
CTRL_list = all_sample_input_table[,2]
TP_list = apply(all_sample_input_table, 1, function(x) unlist(strsplit(x[1], split='_'))[2])
ID_list = apply(all_sample_input_table, 1, function(x) unlist(strsplit(x[1], split='_'))[1])

### get LM scale factors
LM_sf = read.table('s3norm_SF.6TP.txt', header=TRUE)
### 
method_num = dim(LM_sf)[2] / length(IP_list)
sample_num = length(IP_list)

### get LM adjusted signal matrix
signal_mat_colnames = c()
signal_mat = c()
### start
for (i in c(1:sample_num)){
	### get unadjusted signal vector
	signal_file = paste('merge_peak_tab/', IP_list[i], '.allmerged.sorted.signal.fc.tab', sep='')
	signal_vec = scan(signal_file)
	sf_tmp = LM_sf[5, i]
	print(sf_tmp)
	### LM adjustment
	signal_vec = signal_vec * sf_tmp
	#signal_vec = (signal_vec * sf_tmp[2] + sf_tmp[1])
	### cbind matrix
	signal_mat = cbind(signal_mat, signal_vec)
}

### get signal_mat sample IDs
ID_TP = apply(cbind(TP_list, ID_list), 1, function(x) paste(x[1], x[2], sep='_'))
colnames(signal_mat) = ID_TP

### get peak name sorted bed
peak_name_sort_pks = as.vector(read.table('allmerged.peakname_sorted.bed.txt', header=FALSE))
rownames(signal_mat) = peak_name_sort_pks[,1]

### write output matrix
write.table(signal_mat, 'qPCR_cor_6TP/signal_mat.6TP.txt', col.names=T, row.names=T, sep='\t', quote=F)

### Mclust 
cluster_num = 4
fit = Mclust(log2(signal_mat+1), G = cluster_num)#, modelNames = 'VVV')
fit_mean = c()
fit_id_mclust = fit$classification
for (i in c(1:cluster_num)){
	sig_tmp = signal_mat[fit$classification==i,]
	sig_tmp_mean = mean(sig_tmp)
	fit_mean[i] = sig_tmp_mean
}

plot_cluster_rank = order(fit_mean)

k=0
signal_mat_mclust_cluster_newid = c()
signal_mat_mclust_cluster = c()

for (i in plot_cluster_rank){
	k = k+1
	used_id_tmp = fit_id_mclust==i
	pk_num_in_cluster = sum(used_id_tmp)
	signal_mat_mclust_cluster = rbind(signal_mat_mclust_cluster, signal_mat[used_id_tmp,])
	signal_mat_mclust_cluster_newid = c(signal_mat_mclust_cluster_newid, rep(k, pk_num_in_cluster))
}


### plot pheatmap
png('qPCR_cor_6TP/signal_mat_Mclust.6TP.png')
my_colorbar = colorRampPalette(c('white', 'red'))(n = 128)
pheatmap(signal_mat_mclust_cluster, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

pdf('qPCR_cor_6TP/signal_mat_Mclust.6TP.pdf')
my_colorbar = colorRampPalette(c('white', 'red'))(n = 128)
pheatmap(signal_mat_mclust_cluster, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()



