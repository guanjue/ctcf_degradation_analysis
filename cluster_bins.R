cp IDEAS_run_linear/run_IDEAS_result/run_IDEAS.wg.sort.state.bed ./

tail -n+2 run_IDEAS.wg.sort.state.bed | awk -F '\t' -v OFS='\t' '{if ($4!="0_0_0_0") print $0}' | sort -k1,1 -k2,2n > run_IDEAS.wg.sort.state.pk.bed


bedtools intersect -a 201803106_AllCtCFPeaks_DiffBind.bed -b run_IDEAS.wg.sort.state.pk.bed -wa -u > diffbind.ideas.pk.bed

bedtools merge -i  run_IDEAS.wg.sort.state.pk.bed >  run_IDEAS.wg.sort.state.pk.merge.bed



.2r_nbp.fisher_p.200bp.txt


paste mm9_200.sort.bed $ct'.2r_nbp.fisher_p.200bp.txt' > $ct'.2r_nbp.fisher_p.200bp.bedgraph'

time /storage/home/gzx103/group/software/ucsc/bedGraphToBigWig $ct'.2r_nbp.fisher_p.200bp.normed.txt' /storage/home/gzx103/group/projects/vision/input_norm/mm9.chrom.sizes $ct'.S3norm.bw'

paste WT.2r_nbp.fisher_p.200bp.normed.txt 0hr.2r_nbp.fisher_p.200bp.normed.txt 4hr.2r_nbp.fisher_p.200bp.normed.txt 6hr.2r_nbp.fisher_p.200bp.normed.txt > mm9.200bp.signal.mat.txt

paste mm9_200.sort.bed signal_mat_index_vec.txt mm9.200bp.signal.mat.txt > mm9_200.sort.indexset.sig.txt

cat mm9_200.sort.indexset.sig.txt | awk -F '\t' -v OFS='\t' '{if ($4=="6_6_6_6") print $0}' > mm9_200.sort.indexset.sig.high_stable.txt

bedtools merge -i mm9_200.sort.indexset.sig.high_stable.txt > mm9_200.sort.indexset.sig.high_stable.merge.bed


signal_mat = read.table('mm9.200bp.signal.mat.txt', header=F)
signal_mat = as.matrix(signal_mat)
signal_mat_log2 = log2(signal_mat+0.1)


png('hist.S3norm.log2.png')
hist(signal_mat_log2, breaks = 50)
dev.off()

signal_mat_log2_max = apply(signal_mat_log2, 1, max)
png('hist.S3norm.log2.max.png')
hist(signal_mat_log2_max, breaks = 50)
dev.off()

signal_mat_log2_mean = apply(signal_mat_log2, 1, mean)
png('hist.S3norm.log2.mean.png')
hist(signal_mat_log2_mean, breaks = 50)
dev.off()


d1 = scan('WT.rep1.nbp_1r_bgadj.txt')
d2 = scan('WT.rep2.nbp_1r_bgadj.txt')

set.seed(2018)
sample_id = sample(length(d1)[1], 500000)

d1_s = d1[sample_id]
d2_s = d2[sample_id]

png('check.scatter.png')
plot(log2(d1_s+0.1), log2(d2_s+0.1))
dev.off()


library(mclust)

png('hist.all.png')
par(mfrow=c(2,2))
hist(signal_mat_log2[,1], breaks = 50)
hist(signal_mat_log2[,2], breaks = 50)
hist(signal_mat_log2[,3], breaks = 50)
hist(signal_mat_log2[,4], breaks = 50)
dev.off()

#noise_thresh = 2
#set.seed(2018)
signal_mat_vec = as.vector(signal_mat)
signal_mat_p_vec = 10^(-signal_mat_vec)
signal_mat_p_fdr_vec = p.adjust(signal_mat_p_vec, 'fdr')
signal_mat_vec_pk = signal_mat_vec[signal_mat_p_fdr_vec<0.01]

fdr_sig_thresh = log2(min(signal_mat_vec_pk))

#sample_id = sample(dim(signal_mat_log2)[1], 300000)
signal_mat_log2_vec = as.vector(signal_mat_log2)

signal_mat_log2_vec_high = signal_mat_log2_vec[signal_mat_log2_vec>=fdr_sig_thresh]

print(length(signal_mat_log2_vec_high))


### 2nd GMM
set.seed(2018)
mod_all <- densityMclust(signal_mat_log2_vec_high)

summary(mod_all)
attributes(mod_all)

cluster_id = mod_all$classification
cluster_mean = mod_all$parameters$mean
print('2nd GMM cluster means: ')
print(cluster_mean)

rainbow_cp = rev(rainbow(length(cluster_mean)))

png('gmm.density.png')
par(mfrow=c(2,2))
plot(mod_all, what = "density", data = signal_mat_log2_vec_high, breaks = 50)
#plotDensityMclust1(mod_all, data = signal_mat_log2_vec_high, hist.col = "lightgrey", hist.border = "white",  breaks = "Sturges", type = "persp")
#plot(mod_all, what = "density", data = signal_mat_log2_vec_high, drawlabels = FALSE, points.pch = 20)
for (i in c(1:length(cluster_mean))){
	print(i)
	x_input = seq(2,9, 0.1)
	cp_i_mean = mod_all$parameters$mean[i]
	cp_i_sd = (mod_all$parameters$variance$sigmasq[i])^0.5
	cp_i_pro = mod_all$parameters$pro[i]
	lines(x_input, cp_i_pro * dnorm(x_input, mean=cp_i_mean, sd=cp_i_sd), col=rainbow_cp[i])
}
plot(mod_all, what = "BIC")
#plot(mod_all, what = "diagnostic", type = "cdf")
#plot(mod_all, what = "diagnostic", type = "qq")
dev.off()



### get each cluster threshold
gmm_2nd_thresh = c()
gmm_2nd_thresh[1] = fdr_sig_thresh
for (i in c(2:(length(cluster_mean)))){
	print(i)
	gmm_2nd_thresh[i] = min(signal_mat_log2_vec_high[cluster_id==i])
}
print(gmm_2nd_thresh)

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

write.table(index_count, 'index_count_all.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
write.table(index_count_enrich, 'index_count_enrich.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')












signal_mat_log2_s_h = signal_mat_log2_s[cluster_id==4]
mod_pk <- densityMclust(signal_mat_log2_s_h)
png('gmm.density.pk.png')
par(mfrow=c(2,2))
plot(mod_pk, what = "density", data = signal_mat_log2_s_h, breaks = 50)
plot(mod_pk, what = "BIC")
#plot(mod_pk, what = "diagnostic", type = "cdf")
#plot(mod_pk, what = "diagnostic", type = "qq")
dev.off()

summary(mod_pk)






library(mixtools)
gmm_k = 5
mixmdl = normalmixEM(signal_mat_log2_s, k = gmm_k)
png('tp.sig.log2.gmm.png')
plot(mixmdl, which = 2)
dev.off()

mixmdl$mu


multmixmodel.sel(signal_mat_log2_s)

signal_mat_log2_s = as.matrix(signal_mat_log2_s)

set.seed(2018)
library(mixtools)
gmm_k = 4
mixmdl = normalmixEM(signal_mat_log2_s[signal_mat_log2_s>2], k = gmm_k)
png('tp.sig.log2.gmm.png')
plot(mixmdl, which = 2)
lines(density(signal_mat_log2_s[signal_mat_log2_s>2]), lty = 2, lwd = 2)
dev.off()


c = rowSums(signal_mat)

signal_mat_label = signal_mat

signal_mat_label = signal_mat_label[c]
signal_mat_label[signal_mat<=2]=0
signal_mat_label[(signal_mat>2)*(signal_mat<=2.6)==1]=0
signal_mat_label[(signal_mat>2.6)*(signal_mat<=3.8)==1]=1
signal_mat_label[(signal_mat>3.8)*(signal_mat<=6.3)==1]=2
signal_mat_label[(signal_mat>6.3)]=3

signal_mat_label_fast = signal_mat_label[sample_id,]
signal_mat_label_vec = apply(signal_mat_label_fast,1,function(x) paste(x[1], x[2], x[3], x[4], sep='_'))

tl = as.matrix(table(signal_mat_label_vec))

tl[order(rownames(tl))]
