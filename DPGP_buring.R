label=read.table('output_all_1001_percentage_optimal_clustering.txt', header=T)
sigmat0 = read.table('rep_sig_TSnorm_mat_ave.percent.r.DPGP.sub.txt', header=T)

#label=read.table('output_all_1001_percentage_optimal_clustering.txt', header=T)
#sigmat0 = read.table('rep_sig_s3norm_mat_ave.percent.r.DPGP.sub.txt', header=T)

label_vec = label[,1]
dr_nodrb = sigmat0[,-c(1,2)]


library(MASS)
library(pheatmap)
set.seed(2018)



label_vec_factor = factor(label_vec)
qda_input_mat = as.data.frame(dr_nodrb)
qda_input_mat = cbind(qda_input_mat, label_vec_factor)
snapshot_qda_class_mat = c()


rare_clust = rownames(table(label_vec_factor))[table(label_vec_factor)<50]
for (c in rare_clust){
label_vec_factor = as.matrix(label_vec_factor)
label_vec_factor[label_vec_factor==c] = 'X_X'
label_vec_factor = factor(label_vec_factor)
}


for (i in c(1:100)){
print(cbind(table(label_vec_factor)))
snapshot_qda.fit = qda(label_vec_factor ~ X4A+X6A+X12A+X18A+X24A, data = qda_input_mat)
snapshot_qda.class = predict(snapshot_qda.fit, qda_input_mat)$class
label_vec_factor = snapshot_qda.class

if (i%%1==0){
snapshot_qda_class_mat = cbind(snapshot_qda_class_mat, table(snapshot_qda.class))
print(table(snapshot_qda.class))
rare_clust = rownames(table(label_vec_factor))[table(label_vec_factor)<50]
for (c in rare_clust){
label_vec_factor = as.matrix(label_vec_factor)
label_vec_factor[label_vec_factor==c] = 'X_X'
label_vec_factor = factor(label_vec_factor)
}
}

qda_input_mat = as.data.frame(dr_nodrb)
qda_input_mat = cbind(qda_input_mat, label_vec_factor)
}


sigmat0_all = read.table('rep_sig_TSnorm_mat_ave.percent.all.DPGP.txt', header=T)
ave = as.matrix(read.table('rep_sig_TSnorm_mat_ave.txt', header=TRUE))

qda_input_mat_all = as.data.frame(sigmat0_all[,-c(1,2)])
qda_input_mat_all_class = predict(snapshot_qda.fit, qda_input_mat_all)$class



dr_mclust_rep = log2(ave[order(qda_input_mat_all_class),-c(7:10)]+1)

dr_mclust_plot = dr_mclust_rep
nr=10
dr_mclust_plot[dr_mclust_plot>7]=7
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.DPGP.all.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

for (i in rownames(table(qda_input_mat_all_class))){
	sigmat = log2(ave[qda_input_mat_all_class==i,-c(7:10)]+1)
	print(dim(sigmat))
	pdf(paste('mclust.dpgp.box.', toString(i), '.pdf', sep=''))
	boxplot(as.matrix(sigmat), ylim=c(min(log2(ave[,-1]+1)), max(log2(ave[,-1]+1))) ,outline=FALSE)
	for (j in c(1:500)){
		lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	lines(c(1:(dim(sigmat)[2])), colMeans(sigmat), col='black')
	dev.off()
}




dr_mclust_rep = log2(sigmat0[order(label_vec_factor),-1])

dr_mclust_plot = dr_mclust_rep
nr=10
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.DPGP.all.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

for (i in rownames(table(label_vec_factor))){
	sigmat = log2(sigmat0[label_vec_factor==i,-1])
	print(dim(sigmat))
	pdf(paste('mclust.box.', toString(i), '.pdf', sep=''))
	boxplot(as.matrix(sigmat), ylim=c(min(log2(sigmat0[,-1])), max(log2(sigmat0[,-1]))) ,outline=FALSE)
	for (j in c(1:dim(sigmat)[1])){
		lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	lines(c(1:(dim(sigmat)[2])), colMeans(sigmat), col='black')
	dev.off()
}







