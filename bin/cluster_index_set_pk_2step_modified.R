
signal_mat = read.table('201803106_AllCtCFPeaks_DiffBind.S3norm.nochrUn_random.txt', header=F)
signal_mat = as.matrix(signal_mat)

signal_mat_log2 = log2(signal_mat+0.1)
signal_mat_log2_od = signal_mat_log2

signal_mat_log2_max = apply(signal_mat_log2, 1, max)
png('hist.S3norm.log2.max.png')
hist(signal_mat_log2_max, breaks = 50)
dev.off()

signal_mat_log2_mean = apply(signal_mat_log2, 1, mean)
png('hist.S3norm.log2.mean.png')
hist(signal_mat_log2_mean, breaks = 50)
dev.off()


png('hist.S3norm.log2.0hr.png')
hist(signal_mat_log2[,c(2)], breaks = 50)
dev.off()

#signal_mat_log2 = signal_mat_log2-rowMeans(signal_mat_log2)
signal_mat_log2_fc = signal_mat_log2-(signal_mat_log2[,2])


png('hist.S3norm.log2.png')
hist(signal_mat_log2_fc[,c(1,3,4)], breaks = 50)
dev.off()

#d1 = scan('WT.rep1.nbp_1r_bgadj.txt')
#d2 = scan('WT.rep2.nbp_1r_bgadj.txt')

set.seed(2018)
#sample_id = sample(length(d1)[1], 500000)

#d1_s = d1[sample_id]
#d2_s = d2[sample_id]

#png('check.scatter.png')
#plot(log2(d1_s+0.1), log2(d2_s+0.1))
#dev.off()


library(mclust)

png('hist.all.png')
par(mfrow=c(2,2))
hist(signal_mat_log2_fc[,1], breaks = 50)
hist(signal_mat_log2_fc[,2], breaks = 50)
hist(signal_mat_log2_fc[,3], breaks = 50)
hist(signal_mat_log2_fc[,4], breaks = 50)
dev.off()

#noise_thresh = 2
#set.seed(2018)
#signal_mat_vec = as.vector(signal_mat)
#signal_mat_p_vec = 10^(-signal_mat_vec)
#signal_mat_p_fdr_vec = p.adjust(signal_mat_p_vec, 'fdr')
#signal_mat_vec_pk = signal_mat_vec[signal_mat_p_fdr_vec<0.05]

#fdr_sig_thresh = log2(min(signal_mat_vec))

#sample_id = sample(dim(signal_mat_log2)[1], 300000)

### 1st step GMM cluster 0hr signal
signal_mat_0hr_log2_vec = as.vector(signal_mat_log2[,c(2)])
set.seed(2018)
mod_all_0hr_bic <- densityMclust(signal_mat_0hr_log2_vec)
set.seed(2018)
mod_all_0hr <- densityMclust(signal_mat_0hr_log2_vec, G=3)

cluster_id_0hr = mod_all_0hr$classification
cluster_mean_0hr = mod_all_0hr$parameters$mean
print('2nd GMM cluster means: ')
print(cluster_mean_0hr)

rainbow_cp_0hr = rev(rainbow(length(cluster_mean_0hr)))

png('gmm.density.0hr.png')
par(mfrow=c(2,2))
plot(mod_all_0hr, what = "density", data = signal_mat_0hr_log2_vec, breaks = 50)
#plotDensityMclust1(mod_all, data = signal_mat_log2_vec_high, hist.col = "lightgrey", hist.border = "white",  breaks = "Sturges", type = "persp")
for (i in c(1:length(cluster_mean_0hr))){
        print(i)
        x_input = seq(-4,9, 0.1)
        cp_i_mean = mod_all_0hr$parameters$mean[i]
        cp_i_sd = (mod_all_0hr$parameters$variance$sigmasq[i])^0.5
        cp_i_pro = mod_all_0hr$parameters$pro[i]
        lines(x_input, cp_i_pro * dnorm(x_input, mean=cp_i_mean, sd=cp_i_sd), col=rainbow_cp_0hr[i])
}
plot(mod_all_0hr_bic, what = "BIC")
#plot(mod_all, what = "diagnostic", type = "cdf")
#plot(mod_all, what = "diagnostic", type = "qq")
dev.off()


### get each cluster threshold
signal_mat_index_0hr = signal_mat_log2[,2]
signal_0hr_thresh = min(signal_mat_index_0hr)

gmm_0hr_thresh = c()
gmm_0hr_thresh[1] = signal_0hr_thresh
for (i in c(2:(length(cluster_mean_0hr)))){
        print(i)
        gmm_0hr_thresh[i] = min(signal_mat_0hr_log2_vec[cluster_id_0hr==i])
}
print(gmm_0hr_thresh)


### get background index: '0'
for (i in c(1:(length(cluster_mean_0hr)-1))){
        print(i)
        ### get range id
        used_id_tmp = ( (signal_mat_log2[,2]>=gmm_0hr_thresh[i]) * (signal_mat_log2[,2]<gmm_0hr_thresh[i+1]) ) >0
        signal_mat_index_0hr[used_id_tmp] = i
}

signal_mat_index_0hr[signal_mat_log2[,2]>=gmm_0hr_thresh[length(gmm_0hr_thresh)]] = (length(gmm_0hr_thresh))


############################
############################
############################
### 2nd step GMM
signal_mat_log2_fc_vec = as.vector(signal_mat_log2_fc[,c(1,3,4)])
signal_fc_thresh = min(signal_mat_log2_fc_vec)

set.seed(2018)
mod_all_bic <- densityMclust(signal_mat_log2_fc_vec)
set.seed(2018)
mod_all <- densityMclust(signal_mat_log2_fc_vec, G=3)

summary(mod_all)
attributes(mod_all)

cluster_id = mod_all$classification
cluster_mean = mod_all$parameters$mean
print('2nd GMM cluster means: ')
print(cluster_mean)

rainbow_cp = rev(rainbow(length(cluster_mean)))

png('gmm.density.png')
par(mfrow=c(2,2))
plot(mod_all, what = "density", data = signal_mat_log2_fc_vec, breaks = 50)
#plotDensityMclust1(mod_all, data = signal_mat_log2_vec_high, hist.col = "lightgrey", hist.border = "white",  breaks = "Sturges", type = "persp")
for (i in c(1:length(cluster_mean))){
	print(i)
	x_input = seq(-4,9, 0.1)
	cp_i_mean = mod_all$parameters$mean[i]
	cp_i_sd = (mod_all$parameters$variance$sigmasq[i])^0.5
	cp_i_pro = mod_all$parameters$pro[i]
	lines(x_input, cp_i_pro * dnorm(x_input, mean=cp_i_mean, sd=cp_i_sd), col=rainbow_cp[i])
}
plot(mod_all_bic, what = "BIC")
#plot(mod_all, what = "diagnostic", type = "cdf")
#plot(mod_all, what = "diagnostic", type = "qq")
dev.off()



### get each cluster threshold
gmm_2nd_thresh = c()
gmm_2nd_thresh[1] = signal_fc_thresh
for (i in c(2:(length(cluster_mean)))){
	print(i)
	gmm_2nd_thresh[i] = min(signal_mat_log2_fc_vec[cluster_id==i])
}
print(gmm_2nd_thresh)

signal_mat_index = signal_mat_log2_fc
### get background index: '0'
for (i in c(1:(length(cluster_mean)-1))){
	print(i)
	### get range id
	used_id_tmp = ( (signal_mat_log2_fc>=gmm_2nd_thresh[i]) * (signal_mat_log2_fc<gmm_2nd_thresh[i+1]) ) >0
	signal_mat_index[used_id_tmp] = i
}


### get top pk index
signal_mat_index[signal_mat_log2_fc>=gmm_2nd_thresh[length(gmm_2nd_thresh)]] = length(gmm_2nd_thresh)

signal_mat_index = cbind(signal_mat_index_0hr, signal_mat_index)

### get index
signal_mat_index_vec = apply(signal_mat_index, 1, function(x) paste(toString(x[1]), toString(x[2]), toString(x[4]), toString(x[5]), sep="_"))

write.table(signal_mat_index_vec, 'signal_mat_index_vec.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

### get index count
index_count = as.matrix(table(signal_mat_index_vec))
png('index_set_count_hist.png')
hist(log10(index_count))
dev.off()
index_count_enrich = as.matrix(index_count[index_count[,1]>200,])

write.table(index_count, 'index_count_all.txt', quote=FALSE, col.names=FALSE, row.names=TRUE, sep='\t')
write.table(index_count_enrich, 'index_count_enrich.txt', quote=FALSE, col.names=FALSE, row.names=TRUE, sep='\t')






