
### clustering

library(pheatmap)

library(mclust)


sample_list_file = args[1]
sample_list_file = 'id_list_all.txt'
sample_list = read.table(sample_list_file, header=FALSE)
sample_list = sample_list[-c(9,14),]

tp_type = c()
i=1
for (n in sample_list[,3]){
	print(n)
	if (! n %in% tp_type){
		tp_type[i] = n
		i = i+1
	}
}


rep = as.matrix(read.table('rep_sig_s3norm_mat.txt', header=TRUE))
ave = as.matrix(read.table('rep_sig_s3norm_mat_ave.txt', header=TRUE))

#rep = as.matrix(read.table('rep_sig_s3norm_mat.r.txt', header=TRUE))
#ave = as.matrix(read.table('rep_sig_s3norm_mat_ave.r.txt', header=TRUE))

rep_nopk = as.matrix(read.table('rep_sig_s3norm_mat_nopk.txt', header=TRUE))
ave_nopk = as.matrix(read.table('rep_sig_s3norm_mat_nopk_ave.txt', header=TRUE))

rep = as.matrix(read.table('rep_sig_TSnorm_mat.txt', header=TRUE))
ave = as.matrix(read.table('rep_sig_TSnorm_mat_ave.txt', header=TRUE))

rep_nopk = as.matrix(read.table('rep_sig_TSnorm_mat_nopk.txt', header=TRUE))
ave_nopk = as.matrix(read.table('rep_sig_TSnorm_mat_nopk_ave.txt', header=TRUE))



set.seed(2018)
used_id_nopk = sample(dim(ave_nopk)[1], 5000)


scale_fc = function(x){
	xs = (log2(x))
	return((xs)-(xs[1]))
}


scale_log2 = function(x){
	xs = (log2(x))
	return((xs))
}


ave_fc = t(apply(ave+0.1, 1, scale_fc))
#ave_fc = t(apply(ave+0.1, 1, scale_log2))


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
nr = 10
dr_nodrb = as.matrix(ave_fc[,-c(7:10)])


dr = as.matrix(ave_fc[,])


pdf('fc_hist.pdf', width=7, height=42)
par(mfrow=c(6,1))
hist(dr_nodrb[,-6], breaks=50, xlim=c(-3,8), ylim=c(0,12000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(dr_nodrb[,1], breaks=50, xlim=c(-3,8), ylim=c(0,12000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(dr_nodrb[,2], breaks=50, xlim=c(-3,8), ylim=c(0,12000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(dr_nodrb[,3], breaks=50, xlim=c(-3,8), ylim=c(0,12000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(dr_nodrb[,4], breaks=50, xlim=c(-3,8), ylim=c(0,12000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(dr_nodrb[,5], breaks=50, xlim=c(-3,8), ylim=c(0,12000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
dev.off()



plot_lim = 6
### cluster without drb
fit = kmeans(dr_nodrb, 4)
dr_kmeans = log2(rep[order(fit$cluster),]+1)
dr_kmeans_plot = rbind(dr_kmeans, rep_nopk[used_id_nopk,])
dr_kmeans_plot[dr_kmeans_plot>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.s3.nodrb.', toString(nr), '.pdf', sep=''))
pheatmap(dr_kmeans_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

dr_kmeans_ave = log2(ave[order(fit$cluster),]+1)
dr_kmeans_ave_plot = rbind(dr_kmeans_ave, ave_nopk[used_id_nopk,])
dr_kmeans_ave_plot[dr_kmeans_ave_plot>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.s3.ave.nodrb.', toString(nr), '.pdf', sep=''))
pheatmap(dr_kmeans_ave_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

png(paste('kmean.s3.ave.nodrb.', toString(nr), '.png', sep=''))
pheatmap(rbind(dr_kmeans_ave, ave_nopk[used_id_nopk,]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()



library(mclust)

set.seed(2018)

BIC = mclustBIC(dr_nodrb)
png('bic.png')
plot(BIC)
dev.off()

#fit = Mclust(dr_nodrb, x = BIC)
fit = Mclust(dr_nodrb, G=4)
dr_mclust_ave = log2(ave[order(fit$classification),]+1)

#dr_mclust_plot_ave = rbind(dr_mclust_ave, ave_nopk[used_id_nopk,])
dr_mclust_plot_ave = dr_mclust_ave
plot_lim=100
dr_mclust_plot_ave[dr_mclust_plot_ave>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.s3.ave.logsigclust.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot_ave[,-c(7:10)], color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()


for (i in c(1:dim(table(fit$classification))[1])){
	sigmat = log2(ave[fit$classification==i,-c(7:10)]+1)
	print(dim(sigmat)/dim(ave))
	pdf(paste('mclust.TS_snap.box.', toString(i), '.pdf', sep=''))
	boxplot(as.matrix(sigmat), ylim=c(min(log2(ave+1)), max(log2(ave+1))) ,outline=FALSE)
	if (dim(sigmat)[1]>500){line_num = 500} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	lines(c(1:(dim(sigmat)[2])), colMeans(sigmat), col='black')
	dev.off()
}







dr_mclust_ave = log2(ave[order(fit$classification),]+1)

#dr_mclust_plot_ave = rbind(dr_mclust_ave, ave_nopk[used_id_nopk,])
dr_mclust_plot_ave = dr_mclust_ave

dr_mclust_plot_ave[dr_mclust_plot_ave>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.s3.ave.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot_ave[,-c(7:10)], color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

dr_mclust_rep = log2(rep[order(fit$classification),]+1)

#dr_mclust_plot = rbind(dr_mclust_rep, rep_nopk[used_id_nopk,])
dr_mclust_plot = dr_mclust_rep

dr_mclust_plot[dr_mclust_plot>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.s3.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()






write.table(fit$classification, 'mclust.classification.txt', quote=F, sep='\t', col.names=F, row.names=F)

for (i in c(1:dim(table(fit$classification))[1])){
	sigmat = log2(ave[fit$classification==i,]+1)
	print(dim(sigmat))
	pdf(paste('mclust.box.logsig.', toString(i), '.pdf', sep=''))
	boxplot(as.matrix(sigmat), ylim=c(min(log2(ave+1)), max(log2(ave+1))) ,outline=FALSE)
	lines(c(1:(dim(sigmat)[2]-4)), colMeans(sigmat)[-c(7:10)], col='black')
	lines(c(7,8), colMeans(sigmat)[c(7,8)], col='red')
	lines(c(9,10), colMeans(sigmat)[c(9,10)], col='blue')
	dev.off()
}




paste allmerged.sorted.bed.txt rep_sig_s3norm_mat.txt > rep_sig_s3norm_mat.txt.bed
paste allmerged.sorted.bed.txt rep_sig_s3norm_mat_ave.txt > rep_sig_s3norm_mat_ave.txt.bed
paste allmerged.sorted.bed.txt rep_sig_TSnorm_mat.txt > rep_sig_TSnorm_mat.txt.bed
paste allmerged.sorted.bed.txt rep_sig_TSnorm_mat_ave.txt > rep_sig_TSnorm_mat_ave.txt.bed



wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz
gunzip gencode.vM1.annotation.gtf.gz
cat gencode.vM1.annotation.gtf | awk '{if ($3=="gene") print $0}' > gencode.vM1.annotation.gene.gtf
cat gencode.vM1.annotation.gene.gtf | awk -F '"' '{if ($6=="protein_coding") print $0}' > gencode.vM1.annotation.gene.pc.gtf
cat gencode.vM1.annotation.gene.pc.gtf | awk '{print $1,$4,$5,$2,$7,$9}' > gencode.vM1.annotation.gene.pc.gtf.bed
sort -k1,1 -k2,2n gencode.vM1.annotation.gene.pc.gtf.bed > gencode.vM1.annotation.gene.pc.gtf.sort.bed
cat gencode.vM1.annotation.gene.gtf | awk '{print $1,$4,$5,$2,$7,$9}' > gencode.vM1.annotation.gene.gtf.bed
sort -k1,1 -k2,2n gencode.vM1.annotation.gene.gtf.bed > gencode.vM1.annotation.gene.gtf.sort.bed

for i in {1..9}
do
	echo $i
	paste allmerged.sorted.bed mclust.classification.txt | awk -F '\t' -v OFS='\t' -v cluster=$i '{if ($5==cluster) print $0}' > 'allmerged.sorted.c'$i'.bed'
	sort -k1,1 -k2,2n 'allmerged.sorted.c'$i'.bed' > 'allmerged.sorted.c'$i'.sort.bed'
	bedtools intersect -a 'allmerged.sorted.c'$i'.sort.bed' -b gencode.vM1.annotation.gene.pc.gtf.sort.bed -wa -u > 'allmerged.sorted.c'$i'.sort.pc_gene.bed'
	a=$(wc -l 'allmerged.sorted.c'$i'.sort.pc_gene.bed' | awk -F ' ' '{print $1}' )
	b=$(wc -l 'allmerged.sorted.c'$i'.sort.bed' | awk -F ' ' '{print $1}')
	echo $a
	echo $b
	echo echo $(( b / a )) | sed 's/..$/.&/'
done








log2_sig = log2(ave[,-c(7:10)]+0.1)
pdf('signal_hist.pdf', width=7, height=49)
par(mfrow=c(7,1))
hist(log2_sig, breaks=50, xlim=c(-3,7), ylim=c(0,20000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(log2_sig[,1], breaks=50, xlim=c(-3,7), ylim=c(0,20000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(log2_sig[,2], breaks=50, xlim=c(-3,7), ylim=c(0,20000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(log2_sig[,3], breaks=50, xlim=c(-3,7), ylim=c(0,20000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(log2_sig[,4], breaks=50, xlim=c(-3,7), ylim=c(0,20000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(log2_sig[,5], breaks=50, xlim=c(-3,7), ylim=c(0,20000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
hist(log2_sig[,6], breaks=50, xlim=c(-3,7), ylim=c(0,20000))
abline(v=0, col='red', lwd=1.5, lty=2)
box()
dev.off()

set.seed(2018)
not_used_col = -6
mod_all_bic <- densityMclust(as.vector(log2_sig[,1]), model='V')
set.seed(2018)
mod_all <- densityMclust(as.vector(log2_sig[,1]), G=3, model='V')

cluster_id = mod_all$classification
cluster_mean = mod_all$parameters$mean
print('2nd GMM cluster means: ')
print(cluster_mean)

rainbow_cp = rev(rainbow(length(cluster_mean)))

rainbow_cp = c('blue', 'green', 'red')

pdf('gmm.density.pdf', width=14, height=7)
par(mfrow=c(1,2))
plot(mod_all, what = "density", data = as.vector(log2_sig[,1]), breaks = 50)
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

set.seed(2018)

gmm_2nd_thresh = c()
f = function(x) dnorm(x, m=cluster_mean[1], sd=mod_all$parameters$variance$sigmasq[1]) * mod_all$parameters$pro[1] - dnorm(x, m=cluster_mean[2], sd=mod_all$parameters$variance$sigmasq[2]) * mod_all$parameters$pro[2]
gmm_2nd_thresh[1] = uniroot(f, interval=c(3, 5))$root
f = function(x) dnorm(x, m=cluster_mean[3], sd=mod_all$parameters$variance$sigmasq[3]) * mod_all$parameters$pro[3] - dnorm(x, m=cluster_mean[2], sd=mod_all$parameters$variance$sigmasq[2]) * mod_all$parameters$pro[2]
gmm_2nd_thresh[2] = uniroot(f, interval=c(4, 6))$root

### signal three clusters
cluster1 = as.vector(log2_sig[,1]) <= gmm_2nd_thresh[1]
cluster2 = (as.vector(log2_sig[,1]) > gmm_2nd_thresh[1]) & (as.vector(log2_sig[,1]) <= gmm_2nd_thresh[2])
cluster3 = as.vector(log2_sig[,1]) > gmm_2nd_thresh[2]

cluster123 = cbind(cluster1,cluster2,cluster3)


library(mclust)

plot_lim = 7
cluster_num_vec = c(1,5,5)
label_vec = rep('0', dim(cluster123)[1])

for (i in c(1:3)){
set.seed(2018)

cluster_ci = cluster123[,i]
dr_nodrb_ci = dr_nodrb[cluster_ci,]
BIC = mclustBIC(dr_nodrb_ci, modelNames=c('EII', 'VII', 'EEI', 'VVV'))
png(paste('bic', toString(i), '.png', sep=''))
plot(BIC)
dev.off()

fit = Mclust(dr_nodrb_ci, x = BIC)
fit = Mclust(dr_nodrb_ci, G=cluster_num_vec[i])
print(summary(fit))

dr_mclust_ave = log2(ave[cluster_ci,][order(fit$classification),]+1)

dr_mclust_plot_ave = rbind(dr_mclust_ave, ave_nopk[used_id_nopk,])

dr_mclust_plot_ave[dr_mclust_plot_ave>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.s3.ave.', toString(nr), '.', toString(i), '.pdf', sep=''))
pheatmap(dr_mclust_plot_ave, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

dr_mclust_rep = log2(rep[cluster_ci,][order(fit$classification),]+1)

dr_mclust_plot = rbind(dr_mclust_rep, rep_nopk[used_id_nopk,])

dr_mclust_plot[dr_mclust_plot>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.s3.', toString(nr), '.', toString(i), '.pdf', sep=''))
pheatmap(dr_mclust_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

label_tmp = apply(as.matrix(fit$classification), 1, function(x) paste(toString(i), x, sep='_'))
label_vec[cluster_ci] = label_tmp
}



library(MASS)


label_vec_factor = factor(label_vec)
qda_input_mat = as.data.frame(dr_nodrb)
qda_input_mat = cbind(qda_input_mat, label_vec_factor)
snapshot_qda_class_mat = c()


rare_clust = rownames(table(label_vec_factor))[table(label_vec_factor)<500]
for (c in rare_clust){
label_vec_factor = as.matrix(label_vec_factor)
label_vec_factor[label_vec_factor==c] = 'X_X'
label_vec_factor = factor(label_vec_factor)
}

qda_input_mat = as.data.frame(dr_nodrb)
qda_input_mat = cbind(qda_input_mat, label_vec_factor)


for (i in c(1:200)){
if (sum(label_vec_factor=='X_X')<6){
print('small X_X')
qda_input_mat_tmp = qda_input_mat[label_vec_factor!='X_X',]
label_vec_factor = factor(as.matrix(label_vec_factor[label_vec_factor!='X_X']))
qda_input_mat_tmp$label_vec_factor = label_vec_factor
print(table(qda_input_mat_tmp$label_vec_factor))
snapshot_qda.fit = qda(label_vec_factor ~ X4A+X6A+X12A+X18A+X24A, data = qda_input_mat_tmp)
} else {
snapshot_qda.fit = qda(label_vec_factor ~ X4A+X6A+X12A+X18A+X24A, data = qda_input_mat)
}
snapshot_qda.class = predict(snapshot_qda.fit, qda_input_mat)$class
label_vec_factor = snapshot_qda.class
if (i%%1==0){
#snapshot_qda_class_mat = cbind(snapshot_qda_class_mat, table(snapshot_qda.class))
print(table(snapshot_qda.class))

rare_clust = rownames(table(label_vec_factor))[table(label_vec_factor)<500]
for (c in rare_clust){
label_vec_factor = as.matrix(label_vec_factor)
label_vec_factor[label_vec_factor==c] = 'X_X'
label_vec_factor = factor(label_vec_factor)
}
}
qda_input_mat = as.data.frame(dr_nodrb)
qda_input_mat = cbind(qda_input_mat, label_vec_factor)
}


label_vec_factor0 = label_vec_factor

label_vec_factor = as.matrix(label_vec_factor0)
label_vec_factor[label_vec_factor=='X_X'] = '1_0'

dr_mclust_ave = log2(ave[order(label_vec_factor),]+1)

dr_mclust_plot_ave = rbind(dr_mclust_ave, ave_nopk[used_id_nopk,])
dr_mclust_plot_ave = rbind(dr_mclust_ave)

dr_mclust_plot_ave[dr_mclust_plot_ave>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.TS_snap.ave.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot_ave[,-c(7:10)], color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

dr_mclust_rep = log2(rep[order(label_vec_factor),]+1)

dr_mclust_plot = rbind(dr_mclust_rep, rep_nopk[used_id_nopk,])

dr_mclust_plot[dr_mclust_plot>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.TS_snap.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot[-c(1:table(label_vec_factor)[1]),], color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

dr_mclust_plot[dr_mclust_plot>plot_lim]=plot_lim
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('mclust.TS_snap.all.', toString(nr), '.pdf', sep=''))
pheatmap(dr_mclust_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()



for (i in rownames(table(label_vec_factor))){
	sigmat = log2(ave[label_vec_factor==i,-c(7:10)]+1)
	print(dim(sigmat)/dim(ave))
	pdf(paste('mclust.TS_snap.box.', toString(i), '.pdf', sep=''))
	boxplot(as.matrix(sigmat), ylim=c(min(log2(ave+1)), max(log2(ave+1))) ,outline=FALSE)
	if (dim(sigmat)[1]>500){line_num = 500} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	lines(c(1:(dim(sigmat)[2])), colMeans(sigmat), col='black')
	dev.off()
}





library(ggplot2)
library(ggpubr)
my_comparisons = list( c("X0A", "X4A"), c("X4A", "X6A"), c("X6A", "X12A"), c("X12A", "X18A"), c("X18A", "X24A") )
my_test = list( alternative ='greater', alternative ='greater',alternative ='greater',alternative ='greater',alternative ='greater' )
my_test = list( alternative ='greater' )

for (i in rownames(table(label_vec_factor))){
print(i)
	sigmat = log2(ave[label_vec_factor==i,]+1)
	print(dim(sigmat))
	sig_df = as.data.frame(c())
	for (j in c(1:dim(sigmat)[2])){
		sig_df = rbind(sig_df, cbind(sigmat[,j], rep(colnames(sigmat)[j], dim(sigmat)[1])))
	}
	sig_df = as.data.frame(sig_df)
	colnames(sig_df) = c('sig', 'set')
	sig_df[,1] = as.numeric(as.character(sig_df[, 1]))

	pdf(paste('mclust.box.', toString(i), '.pdf', sep=''))
	abc = ggplot(sig_df, aes(x = set, y = sig), ylim=c(0,7)) + 
	#ggline(sig_df, x = 'set', y = 'sig', add = "mean_se") + 
	geom_boxplot() +
	stat_compare_means(comparisons = my_comparisons, method.args = my_test)
	print(abc)
	dev.off()
}





s3_label_vec_factor = label_vec_factor


TS_label_vec_factor = label_vec_factor


adjustedRandIndex(s3_label_vec_factor, TS_label_vec_factor)

dif_mat = c()
for (c in rownames(table(s3_label_vec_factor))){
tc = cbind(table(TS_label_vec_factor[s3_label_vec_factor==c]))#/sum(s3_label_vec_factor==c))
print(tc)
dif_mat = cbind(dif_mat, tc)
}

colnames(dif_mat) = rownames(table(s3_label_vec_factor))

dif_mat




