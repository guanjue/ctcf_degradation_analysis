### get parameters
args = commandArgs(trailingOnly=TRUE)

sample_list_file = args[1]
sample_list_file = 'id_list_all.txt'
sample_list = read.table(sample_list_file, header=FALSE)

tp_type = c()
i=1
for (n in sample_list[,3]){
	print(n)
	if (! n %in% tp_type){
		tp_type[i] = n
		i = i+1
	}
}

peak_type = read.table('qtPCR_regions.sorted.bed.txt', header=FALSE)[,5]

###########
s3norm = function(ref, tar, pkid, bgid, A, B, converge_thresh){
	ref_pk_m = mean(ref[pkid])
	ref_bg_m = mean(ref[bgid])
	for (i in c(1:100)){
		fb = ref_bg_m * mean(tar[pkid]^B) - ref_pk_m * mean(tar[bgid]^B)
		dfb = ref_bg_m * mean(log(tar[pkid]) * tar[pkid]^B) - ref_pk_m * mean(log(tar[bgid]) * tar[bgid]^B)
		### next step
		B = B - fb/dfb
		A = ref_bg_m / mean(tar[bgid]^B)
		print(paste('Iteration:',toString(i)))
		print(paste('A:', A))
		print(paste('B:', B))

		last_AB = c(A, B)
		### converge
		if (abs(fb/dfb) < converge_thresh){
			print('converged!')
			used_AB = c(A, B)
			break
		}
	}
	if (abs(fb/dfb)>= converge_thresh){
		print('NOT converged...')
		used_AB = last_AB		
	}
	return(used_AB)
}
###########

###########
s3norm_nonpkbg = function(ref, tar, ref_nonpkbg, tar_nonpkbg, A, B, converge_thresh){
	ref_pk_m = mean(ref)
	ref_bg_m = mean(ref_nonpkbg)
	for (i in c(1:100)){
		fb = ref_bg_m * mean(tar^B) - ref_pk_m * mean(tar_nonpkbg^B)
		dfb = ref_bg_m * mean(log(tar) * tar^B) - ref_pk_m * mean(log(tar_nonpkbg) * tar_nonpkbg^B)
		### next step
		B = B - fb/dfb
		A = ref_bg_m / mean(tar_nonpkbg^B)
		print(paste('Iteration:',toString(i)))
		print(paste('A:', A))
		print(paste('B:', B))

		last_AB = c(A, B)
		### converge
		if (abs(fb/dfb) < converge_thresh){
			print('converged!')
			used_AB = c(A, B)
			break
		}
	}
	if (abs(fb/dfb)>= converge_thresh){
		print('NOT converged...')
		used_AB = last_AB		
	}
	return(used_AB)
}
###########

###########
s3norm_nonpkbg_matchmean = function(ref, tar, ref_nonpkbg, tar_nonpkbg){
	sig1_pk_mean = mean(ref)
	sig1_bg_mean = mean(ref_nonpkbg)
	sig2_pk_mean = mean(tar)
	sig2_bg_mean = mean(tar_nonpkbg)
	B = (sig1_pk_mean - sig1_bg_mean) / (sig2_pk_mean - sig2_bg_mean)
	A = 2^(sig1_pk_mean - B * sig2_pk_mean)
	used_AB = c(A, B)
	print('used: ')
	print(used_AB)
	print(sig1_pk_mean)
	print(sig1_bg_mean)
	print(sig2_pk_mean)
	print(sig2_bg_mean)
	return(used_AB)
}
###########


### get initial matrix
d_qtpcr_mat = c()

d_chipseq_mat = c()
d_chipseq_input_mat = c()

d_nonpk_mat = c()
d_nonpk_input_mat = c()

output_vec = c()

tp_vec_all = sample_list[,3]

### get signal matrix
for (i in c(1:length(tp_type))){
	rep_list = sample_list[sample_list[,3] == tp_type[i],1]
	input_list = sample_list[sample_list[,3] == tp_type[i],2][1]
	print(sample_list[,3])
	for (j in c(1:length(rep_list))){
		rep1 = paste(rep_list[j], '.qtPCR_regions.sorted.rc_per_bp.tab', sep='')
		rep1_nonpk = paste(rep_list[j], '.mm9_500.sorted.rc_per_bp.tab', sep='')
		d_chipseq_1 = scan(rep1)*1000
		d_chipseq_mat = cbind(d_chipseq_mat, d_chipseq_1)
		d_chipseq_1_nonpk = scan(rep1_nonpk)*1000
		d_nonpk_mat = cbind(d_nonpk_mat, d_chipseq_1_nonpk)
		input = paste(input_list, '.qtPCR_regions.input_sig.macs.tab', sep='')
		input_nonpk = paste(input_list, '.mm9_500.input_sig.macs.tab', sep='')
		qtpcr = paste('qtPCR_sig.', tp_type[i], '.txt', sep='')
		output = paste(tp_type[i], '.chipseq_vs_qtpcr.png', sep='')
		output_vec[i] = output
		### read input files
		d_qtpcr = scan(qtpcr)*1000
		d_chipseq_input = scan(input)*1000
		d_chipseq_input_nonpk = scan(input_nonpk)*1000
		### get signal matrix
		d_qtpcr_mat = cbind(d_qtpcr_mat, d_qtpcr)
		d_chipseq_input_mat = cbind(d_chipseq_input_mat, d_chipseq_input)
		d_nonpk_input_mat = cbind(d_nonpk_input_mat, d_chipseq_input_nonpk)
	}
}

### plot ChIP-seq replicates signal at qtPCR regions
lim_min_chipseq = min(d_chipseq_mat)
lim_max_chipseq = max(d_chipseq_mat)
color_list = c('purple', 'brown', 'blue', 'red', 'green', 'orange', 'black')

### extract overall scale model for qtPCR
#( d_qtpcr_mat - mean(d_qtpcr_mat) ) / sd(d_qtpcr_mat) ~ ( d_chipseq_i_mat - mean(d_chipseq_i_mat) ) / sd(d_chipseq_i_mat)
d_chipseq_mat_fc = d_chipseq_mat / d_chipseq_input_mat

d_qtpcr_mat_scale = ( d_qtpcr_mat - mean(d_qtpcr_mat) ) / sd(d_qtpcr_mat) * sd(d_chipseq_mat_fc) + mean(d_chipseq_mat_fc) 

### plot qtPCR replicates signal at qtPCR regions
lim_min_chipseq = min(d_chipseq_mat_fc)
lim_max_chipseq = max(d_chipseq_mat_fc)


ref_nonpkbg_0tp_rep1 = scan('/storage/home/g/gzx103/scratch/ctcf/6_time_points/1579_0A.mm9_500.sorted.rc_per_bp.tab')*1000
ref_nonpkbg_0tp_rep2 = scan('/storage/home/g/gzx103/scratch/ctcf/6_time_points/1580_0A.mm9_500.sorted.rc_per_bp.tab')*1000
ref_nonpkbg_0tp_input = scan('/storage/home/g/gzx103/scratch/ctcf/6_time_points/1578_0A.mm9_500.input_sig.macs.tab')*1000


AB_mat = c()
AB_average_mat = c()

AB_TS_mat = c()
AB_TSaverage_mat = c()


### S3norm vs TSnorm and plot
for (i in c(1:length(tp_type))){
	print(i)
	smallnum = 1
	output_i = output_vec[i]

	used_colid = tp_vec_all == tp_type[i]
	used_colid_ref = tp_vec_all == tp_type[1]

	d_qtpcr_mat_scale_i = rowMeans(d_qtpcr_mat_scale[,used_colid])

	### get for ground signals
	d_chipseq_mat_FC_tmp = (d_chipseq_mat[,used_colid]+smallnum)/(d_chipseq_input_mat[,used_colid]+smallnum)
	d_chipseq_average_mat_FC_i = (rowMeans(d_chipseq_mat[,used_colid])+smallnum)/(rowMeans(d_chipseq_input_mat[,used_colid])+smallnum)

	d_chipseq_mat_FC_j_cor_vec = c()
	for (j in c(1:dim(d_chipseq_mat_FC_tmp)[2])){
		d_chipseq_mat_FC_j = d_chipseq_mat_FC_tmp[,j]
		d_chipseq_mat_FC_j_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_mat_FC_j)
		d_chipseq_mat_FC_j_cor_vec[j] = d_chipseq_mat_FC_j_cor
	}
	d_chipseq_average_mat_FC_i_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i)

	### get non_peak region signals
	d_nonpk_mat_FC_tar_tmp = (d_nonpk_mat[,used_colid]+smallnum)/(d_nonpk_input_mat[,used_colid]+smallnum)
	print(head((d_nonpk_mat[,used_colid]+smallnum)))
	print(head((d_nonpk_input_mat[,used_colid]+smallnum)))
	print(head((d_nonpk_mat_FC_tar_tmp)))
	print(dim(d_nonpk_mat_FC_tar_tmp))
	d_nonpk_average_mat_FC_i_tar = (rowMeans(d_nonpk_mat[,used_colid])+smallnum)/(rowMeans(d_nonpk_input_mat[,used_colid])+smallnum)
	### get reference non_peak region signals
	d_nonpk_mat_FC_ref_tmp = (d_nonpk_mat[,used_colid_ref]+smallnum)/(d_nonpk_input_mat[,used_colid_ref]+smallnum)
	#d_nonpk_average_mat_FC_i_ref = (rowMeans(d_nonpk_mat[,used_colid_ref])+smallnum)/(rowMeans(d_nonpk_input_mat[,used_colid_ref])+smallnum)
	d_nonpk_average_mat_FC_i_ref = (rowMeans(cbind(ref_nonpkbg_0tp_rep1+ref_nonpkbg_0tp_rep2))+smallnum)/(ref_nonpkbg_0tp_input+smallnum)

	### get TSnorm factors 
	TSnorm_tmp_vec = c()
	d_chipseq_mat_FC_i_TSnorm_mat = c()
	for (j in c(1:dim(d_chipseq_mat_FC_tmp)[2])){
		TSnorm_tmp = 1.0 / mean(d_chipseq_mat_FC_tmp[,j]) * mean(d_qtpcr_mat_scale_i)
		TSnorm_tmp_vec[j] = TSnorm_tmp
		### TSnorm
		d_chipseq_mat_FC_j_TSnorm = d_chipseq_mat_FC_tmp[,j] * TSnorm_tmp
		d_chipseq_mat_FC_i_TSnorm_mat = cbind(d_chipseq_mat_FC_i_TSnorm_mat, d_chipseq_mat_FC_j_TSnorm)
	}
	TSnorm_average = 1.0 / mean(d_chipseq_average_mat_FC_i) * mean(d_qtpcr_mat_scale_i)
	d_chipseq_average_mat_FC_i_TSnorm = d_chipseq_average_mat_FC_i * TSnorm_average
	d_nonpk_average_mat_FC_i_tar_TSnorm = d_nonpk_average_mat_FC_i_tar * TSnorm_average

	### get s3norm factors 
	s3norm_tmp_vec = c()
	d_chipseq_mat_FC_i_s3norm_mat = c()
	used_pk = peak_type!='1Strong-Stable'
	AB_tmp_mat = c()
	for (j in c(1:dim(d_chipseq_mat_FC_tmp)[2])){
		AB_tmp = s3norm_nonpkbg(d_qtpcr_mat_scale_i[used_pk], d_chipseq_mat_FC_tmp[,j][used_pk], d_nonpk_average_mat_FC_i_ref, d_nonpk_mat_FC_tar_tmp[,j], 1.0, 2.0, 1e-5)
		AB_tmp_mat = rbind(AB_tmp_mat, AB_tmp)
		### s3norm
		d_chipseq_mat_FC_j_s3norm = AB_tmp[1] * d_chipseq_mat_FC_tmp[,j] ^ AB_tmp[2]
		d_chipseq_mat_FC_i_s3norm_mat = cbind(d_chipseq_mat_FC_i_s3norm_mat, d_chipseq_mat_FC_j_s3norm)
	}
	AB_average = s3norm_nonpkbg(d_qtpcr_mat_scale_i[used_pk], d_chipseq_average_mat_FC_i[used_pk], d_nonpk_average_mat_FC_i_ref, d_nonpk_average_mat_FC_i_tar, 1.0, 2.0, 1e-5)
	d_chipseq_average_mat_FC_i_s3norm = AB_average[1] * d_chipseq_average_mat_FC_i ^ AB_average[2]
	d_nonpk_mat_FC_i_tar_s3norm = AB_average[1] * (d_nonpk_average_mat_FC_i_tar ^ AB_average[2])
	### s3norm correlation
	d_chipseq_mat_FC_s3norm_j_cor_vec = c()
	for (j in c(1:dim(d_chipseq_mat_FC_i_s3norm_mat)[2])){
		d_chipseq_mat_FC_s3norm_j = d_chipseq_mat_FC_i_s3norm_mat[,j]
		d_chipseq_mat_FC_s3norm_j_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_mat_FC_s3norm_j)
		d_chipseq_mat_FC_s3norm_j_cor_vec[j] = d_chipseq_mat_FC_s3norm_j_cor
	}
	d_chipseq_average_mat_FC_s3norm_i_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i_s3norm)


	#pdf(output_i, width=14, height=14)
	png(output_i, width=800, height=800)
	par(mfrow=c(2,2))
	main_title = ''
	for (k in c(1:length(d_chipseq_mat_FC_j_cor_vec))){
		main_title = paste(main_title, 'rep:', round(d_chipseq_mat_FC_j_cor_vec[k], 3), sep=' ')
	}
	main_title = paste(main_title, 'ave:', round(d_chipseq_average_mat_FC_i_cor, 3), sep=' ')
	### plot original
	plot(d_nonpk_average_mat_FC_i_ref, d_nonpk_average_mat_FC_i_tar, pch=16, col='gray', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), main = main_title)
	for (j in c(1:dim(d_chipseq_mat_FC_tmp)[2])){
		points(d_qtpcr_mat_scale_i, d_chipseq_mat_FC_tmp[,j], pch=16, col=color_list[j])
	} 
	points(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i, pch=16, col='black')
	points(mean(d_qtpcr_mat_scale_i), mean(d_chipseq_average_mat_FC_i), pch=16, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Strong-Stable']), mean(d_chipseq_average_mat_FC_i[peak_type=='Strong-Stable']), pch=15, cex = 1.5, col='red')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Strong-Dynamic']), mean(d_chipseq_average_mat_FC_i[peak_type=='Strong-Dynamic']), pch=15, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Weak']), mean(d_chipseq_average_mat_FC_i[peak_type=='Weak']), pch=15, cex = 1.5, col='blue')
	points(mean(d_qtpcr_mat_scale_i[peak_type!='Strong-Stable']), mean(d_chipseq_average_mat_FC_i[peak_type!='Strong-Stable']), pch=15, cex = 1.5, col='black')
	abline(0, 1, lwd=1.5, col = 'red')

	### plot TSnorm
	plot(d_nonpk_average_mat_FC_i_ref, d_nonpk_average_mat_FC_i_tar_TSnorm, pch=16, col='gray', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), main = main_title)
	for (j in c(1:dim(d_chipseq_mat_FC_i_TSnorm_mat)[2])){
		points(d_qtpcr_mat_scale_i, d_chipseq_mat_FC_i_TSnorm_mat[,j], pch=16, col=color_list[j])
	} 
	points(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i_TSnorm, pch=16, col='black')
	points(mean(d_qtpcr_mat_scale_i), mean(d_chipseq_average_mat_FC_i_TSnorm), pch=16, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Strong-Stable']), mean(d_chipseq_average_mat_FC_i_TSnorm[peak_type=='Strong-Stable']), pch=15, cex = 1.5, col='red')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Strong-Dynamic']), mean(d_chipseq_average_mat_FC_i_TSnorm[peak_type=='Strong-Dynamic']), pch=15, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Weak']), mean(d_chipseq_average_mat_FC_i_TSnorm[peak_type=='Weak']), pch=15, cex = 1.5, col='blue')
	points(mean(d_qtpcr_mat_scale_i[peak_type!='Strong-Stable']), mean(d_chipseq_average_mat_FC_i_TSnorm[peak_type!='Strong-Stable']), pch=15, cex = 1.5, col='black')
	abline(0, 1, lwd=1.5, col = 'red')


	### plot original again
	plot(d_nonpk_average_mat_FC_i_ref, d_nonpk_average_mat_FC_i_tar, pch=16, col='gray', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), main = main_title)
	for (j in c(1:dim(d_chipseq_mat_FC_tmp)[2])){
		points(d_qtpcr_mat_scale_i, d_chipseq_mat_FC_tmp[,j], pch=16, col=color_list[j])
	} 
	points(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i, pch=16, col='black')
	points(mean(d_qtpcr_mat_scale_i), mean(d_chipseq_average_mat_FC_i), pch=16, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Strong-Stable']), mean(d_chipseq_average_mat_FC_i[peak_type=='Strong-Stable']), pch=15, cex = 1.5, col='red')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Strong-Dynamic']), mean(d_chipseq_average_mat_FC_i[peak_type=='Strong-Dynamic']), pch=15, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Weak']), mean(d_chipseq_average_mat_FC_i[peak_type=='Weak']), pch=15, cex = 1.5, col='blue')
	points(mean(d_qtpcr_mat_scale_i[peak_type!='Strong-Stable']), mean(d_chipseq_average_mat_FC_i[peak_type!='Strong-Stable']), pch=15, cex = 1.5, col='black')
	abline(0, 1, lwd=1.5, col = 'red')


	### s3norm
	main_title_s3 = ''
	for (k in c(1:length(d_chipseq_mat_FC_s3norm_j_cor_vec))){
		main_title_s3 = paste(main_title_s3, 'rep:', round(d_chipseq_mat_FC_s3norm_j_cor_vec[k], 3), sep=' ')
	}
	main_title_s3 = paste(main_title_s3, 'ave:', round(d_chipseq_average_mat_FC_s3norm_i_cor, 3), sep=' ')

	### plot s3norm
	plot(d_nonpk_average_mat_FC_i_ref, d_nonpk_mat_FC_i_tar_s3norm, pch=16, col='gray', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), main = main_title_s3)
	for (j in c(1:dim(d_chipseq_mat_FC_i_s3norm_mat)[2])){
		points(d_qtpcr_mat_scale_i, d_chipseq_mat_FC_i_s3norm_mat[,j], pch=16, col=color_list[j])
	} 
	points(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i_s3norm, pch=16, col='black')
	points(mean(d_qtpcr_mat_scale_i), mean(d_chipseq_average_mat_FC_i_s3norm), pch=16, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Strong-Stable']), mean(d_chipseq_average_mat_FC_i_s3norm[peak_type=='Strong-Stable']), pch=15, cex = 1.5, col='red')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Strong-Dynamic']), mean(d_chipseq_average_mat_FC_i_s3norm[peak_type=='Strong-Dynamic']), pch=15, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[peak_type=='Weak']), mean(d_chipseq_average_mat_FC_i_s3norm[peak_type=='Weak']), pch=15, cex = 1.5, col='blue')
	points(mean(d_qtpcr_mat_scale_i[peak_type!='Strong-Stable']), mean(d_chipseq_average_mat_FC_i_s3norm[peak_type!='Strong-Stable']), pch=15, cex = 1.5, col='black')
	abline(0, 1, lwd=1.5, col = 'red')
	dev.off()

	print(d_chipseq_mat_FC_j_cor_vec)
	print(d_chipseq_mat_FC_s3norm_j_cor_vec)
	AB_mat = rbind(AB_mat, AB_tmp_mat)
	AB_average_mat = rbind(AB_average_mat, AB_average)
	AB_TS_mat = rbind(AB_TS_mat, cbind(TSnorm_tmp_vec, rep(1, length(TSnorm_tmp_vec))))
	AB_TSaverage_mat = rbind(AB_TSaverage_mat, c(TSnorm_average, 1))
}


AB_mat = cbind(as.data.frame(AB_mat), sample_list[,3])
colnames(AB_mat) = c('A', 'B', 'tp')
rownames(AB_mat) = sample_list[,1]
write.table(AB_mat, 's3norm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

AB_average_mat = cbind(as.data.frame(AB_average_mat), tp_type)
colnames(AB_average_mat) = c('A', 'B', 'tp')
rownames(AB_average_mat) = tp_type
write.table(AB_average_mat, 'average_s3norm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

AB_TS_mat = cbind(as.data.frame(AB_TS_mat), sample_list[,3])
colnames(AB_TS_mat) = c('A', 'B', 'tp')
rownames(AB_TS_mat) = sample_list[,1]
write.table(AB_TS_mat, 'TSnorm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

AB_TSaverage_mat = cbind(as.data.frame(AB_TSaverage_mat), tp_type)
colnames(AB_TSaverage_mat) = c('A', 'B', 'tp')
rownames(AB_TSaverage_mat) = tp_type
write.table(AB_TSaverage_mat, 'average_TSnorm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)





### get parameters
args = commandArgs(trailingOnly=TRUE)

sample_list_file = args[1]
sample_list_file = 'id_list_all.txt'
sample_list = read.table(sample_list_file, header=FALSE)

tp_type = c()
i=1
for (n in sample_list[,3]){
	print(n)
	if (! n %in% tp_type){
		tp_type[i] = n
		i = i+1
	}
}

### get pk num and get sample id
file1 = paste(sample_list[1, 1], '.allmerged.sorted.rc_per_bp.tab', sep='')
file1 = scan(file1)
pknum = length(file1)
set.seed(2018)
used_id = sample(pknum, 10000)

### read S3norm AB

AB_used = read.table('s3norm_AB.txt', sep='\t', header=TRUE)
AB_average_used = read.table('average_s3norm_AB.txt', sep='\t', header=TRUE)

AB_TS_used = read.table('TSnorm_AB.txt', sep='\t', header=TRUE)
AB_TS_average_used = read.table('average_TSnorm_AB.txt', sep='\t', header=TRUE)

rep_sig_s3norm_mat = c()
rep_sig_TSnorm_mat = c()
rep_sig_s3norm_mat_nonpk = c()
rep_sig_TSnorm_mat_nonpk = c()

rep_sig_s3norm_mat_ave = c()
rep_sig_TSnorm_mat_ave = c()
rep_sig_s3norm_mat_nonpk_ave = c()
rep_sig_TSnorm_mat_nonpk_ave = c()

smallnum = 1
### normalization
for (i in c(1:length(tp_type))){
	rep_list = sample_list[sample_list[,3] == tp_type[i],1]
	input_list = sample_list[sample_list[,3] == tp_type[i],2][1]
	AB_used_tmp = AB_used[AB_used[,3]==tp_type[i],]
	AB_average_used_tmp = AB_average_used[AB_average_used[,3]==tp_type[i],]
	AB_TS_used_tmp = AB_TS_used[AB_TS_used[,3]==tp_type[i],]
	AB_TS_average_used_tmp = AB_TS_average_used[AB_TS_average_used[,3]==tp_type[i],]


	d_chipseq_1_mat_tmp = c()
	d_chipseq_1_mat_nonpk_tmp = c()

	print(sample_list[,3])
	for (j in c(1:length(rep_list))){
		AB_used_tmp_j = AB_used_tmp[j,]
		AB_TS_used_tmp_j = AB_TS_used_tmp[j,]

		rep1 = paste(rep_list[j], '.allmerged.sorted.rc_per_bp.tab', sep='')
		input = paste(input_list, '.allmerged.input_sig.macs.tab', sep='')
		rep1_nonpk = paste(rep_list[j], '.mm9_500.sorted.rc_per_bp.tab', sep='')
		input_nonpk = paste(input_list, '.mm9_500.input_sig.macs.tab', sep='')

		d_chipseq_1 = scan(rep1)*1000
		d_chipseq_1_input = scan(input)*1000
		d_chipseq_1_nonpk = scan(rep1_nonpk)*1000
		d_chipseq_1_nonpk_input = scan(input_nonpk)*1000

		### get FC
		d_chipseq_1_fc = (d_chipseq_1+smallnum) / (d_chipseq_1_input+smallnum)
		d_chipseq_1_nonpk_fc = (d_chipseq_1_nonpk+smallnum) / (d_chipseq_1_nonpk_input+smallnum)
		### s3norm
		d_chipseq_1_fc_s3norm = as.numeric(AB_used_tmp_j[1]) * d_chipseq_1_fc ^ as.numeric(AB_used_tmp_j[2])
		d_chipseq_1_nonpk_fc_s3norm = as.numeric(AB_used_tmp_j[1]) * d_chipseq_1_nonpk_fc ^ as.numeric(AB_used_tmp_j[2])
		### TSnorm
		d_chipseq_1_fc_TSnorm = as.numeric(AB_TS_used_tmp_j[1]) * d_chipseq_1_fc
		d_chipseq_1_nonpk_fc_TSnorm = as.numeric(AB_TS_used_tmp_j[1]) * d_chipseq_1_nonpk_fc

		### save to matrix
		rep_sig_s3norm_mat = cbind(rep_sig_s3norm_mat, d_chipseq_1_fc_s3norm)
		rep_sig_TSnorm_mat = cbind(rep_sig_TSnorm_mat, d_chipseq_1_fc_TSnorm)
		rep_sig_s3norm_mat_nonpk = cbind(rep_sig_s3norm_mat_nonpk, d_chipseq_1_nonpk_fc_s3norm)
		rep_sig_TSnorm_mat_nonpk = cbind(rep_sig_TSnorm_mat_nonpk, d_chipseq_1_nonpk_fc_TSnorm)
		### save to ave matrix
		d_chipseq_1_mat_tmp = cbind(d_chipseq_1_mat_tmp, d_chipseq_1)
		d_chipseq_1_mat_nonpk_tmp = cbind(d_chipseq_1_mat_nonpk_tmp, d_chipseq_1_nonpk)
	}

	d_chipseq_1_ave_fc_tmp = (rowMeans(d_chipseq_1_mat_tmp)+smallnum)/(d_chipseq_1_input+smallnum)
	d_chipseq_1_nonpk_ave_fc_tmp = (rowMeans(d_chipseq_1_mat_nonpk_tmp)+smallnum)/(d_chipseq_1_nonpk_input+smallnum)
	### ave s3
	rep_sig_s3norm_mat_ave_tmp = as.numeric(AB_used_tmp_j[1]) * d_chipseq_1_ave_fc_tmp ^ as.numeric(AB_used_tmp_j[2])
	rep_sig_s3norm_mat_nonpk_ave_tmp = as.numeric(AB_used_tmp_j[1]) * d_chipseq_1_nonpk_ave_fc_tmp ^ as.numeric(AB_used_tmp_j[2])
	rep_sig_s3norm_mat_ave = cbind(rep_sig_s3norm_mat_ave, rep_sig_s3norm_mat_ave_tmp)
	rep_sig_s3norm_mat_nonpk_ave = cbind(rep_sig_s3norm_mat_nonpk_ave, rep_sig_s3norm_mat_nonpk_ave_tmp)
	### ave TS
	rep_sig_TSnorm_mat_ave_tmp = as.numeric(AB_TS_used_tmp_j[1]) * d_chipseq_1_ave_fc_tmp 
	rep_sig_TSnorm_mat_nonpk_ave_tmp = as.numeric(AB_TS_used_tmp_j[1]) * d_chipseq_1_nonpk_ave_fc_tmp 
	rep_sig_TSnorm_mat_ave = cbind(rep_sig_TSnorm_mat_ave, rep_sig_TSnorm_mat_ave_tmp)
	rep_sig_TSnorm_mat_nonpk_ave = cbind(rep_sig_TSnorm_mat_nonpk_ave, rep_sig_TSnorm_mat_nonpk_ave_tmp)
}

colnames(rep_sig_s3norm_mat) = sample_list[,3]
colnames(rep_sig_TSnorm_mat) = sample_list[,3]
colnames(rep_sig_s3norm_mat_nonpk) = sample_list[,3]
colnames(rep_sig_TSnorm_mat_nonpk) = sample_list[,3]


colnames(rep_sig_s3norm_mat_ave) = tp_type
colnames(rep_sig_TSnorm_mat_ave) = tp_type
colnames(rep_sig_s3norm_mat_nonpk_ave) = tp_type
colnames(rep_sig_TSnorm_mat_nonpk_ave) = tp_type

bed = read.table('allmerged.sorted.bed', header=FALSE)

colnames(bed) = c('chr', 'start','end','pk')

bed_rep_sig_s3norm_mat = cbind(bed, rep_sig_s3norm_mat)
bed_rep_sig_s3norm_mat_ave = cbind(bed, rep_sig_s3norm_mat_ave)
bed_rep_sig_TSnorm_mat = cbind(bed, rep_sig_TSnorm_mat)
bed_rep_sig_TSnorm_mat_ave = cbind(bed, rep_sig_TSnorm_mat_ave)

write.table(rep_sig_s3norm_mat, 'rep_sig_s3norm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_s3norm_mat_nonpk, 'rep_sig_s3norm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_s3norm_mat_ave, 'rep_sig_s3norm_mat_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_s3norm_mat_nonpk_ave, 'rep_sig_s3norm_mat_nopk_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(rep_sig_TSnorm_mat, 'rep_sig_TSnorm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_TSnorm_mat_nonpk, 'rep_sig_TSnorm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_TSnorm_mat_ave, 'rep_sig_TSnorm_mat_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_TSnorm_mat_nonpk_ave, 'rep_sig_TSnorm_mat_nopk_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(bed_rep_sig_s3norm_mat, 'S3norm_mat_rep.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(bed_rep_sig_s3norm_mat_ave, 'S3norm_mat_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(bed_rep_sig_TSnorm_mat, 'TSnorm_mat_rep.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(bed_rep_sig_TSnorm_mat_ave, 'TSnorm_mat_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)



library(pheatmap)


sample_list_file = args[1]
sample_list_file = 'id_list_all.txt'
sample_list = read.table(sample_list_file, header=FALSE)

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

rep_nopk = as.matrix(read.table('rep_sig_s3norm_mat_nopk.txt', header=TRUE))
ave_nopk = as.matrix(read.table('rep_sig_s3norm_mat_nopk_ave.txt', header=TRUE))


set.seed(2018)
used_id_nopk = sample(dim(ave_nopk)[1], 5000)


scale_fc = function(x){
	xs = (log2(x))
	return((xs)-(xs[6]))
}

ave_fc = t(apply(ave, 1, scale_fc))


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
dr = as.matrix(ave_fc)
fit = kmeans(dr, 10)


dr_kmeans = log2(rep[order(fit$cluster),]+1)
dr_kmeans_plot = rbind(dr_kmeans, rep_nopk[used_id_nopk,])
dr_kmeans_plot[dr_kmeans_plot>7]=7
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.s3.', toString(nr), '.pdf', sep=''))
pheatmap(dr_kmeans_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

dr_kmeans_ave = log2(ave[order(fit$cluster),]+1)
dr_kmeans_ave_plot = rbind(dr_kmeans_ave, ave_nopk[used_id_nopk,])
dr_kmeans_ave_plot[dr_kmeans_ave_plot>7]=7
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.s3.ave.', toString(nr), '.pdf', sep=''))
pheatmap(dr_kmeans_ave_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

png(paste('kmean.s3.ave.', toString(nr), '.png', sep=''))
pheatmap(rbind(dr_kmeans_ave, ave_nopk[used_id_nopk,]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()





rep = as.matrix(read.table('rep_sig_TSnorm_mat.txt', header=TRUE))
ave = as.matrix(read.table('rep_sig_TSnorm_mat_ave.txt', header=TRUE))

rep_nopk = as.matrix(read.table('rep_sig_TSnorm_mat_nopk.txt', header=TRUE))
ave_nopk = as.matrix(read.table('rep_sig_TSnorm_mat_nopk_ave.txt', header=TRUE))


set.seed(2018)
used_id_nopk = sample(dim(ave_nopk)[1], 5000)


scale_fc = function(x){
	xs = (log2(x))
	return((xs)-(xs[6]))
}

ave_fc = t(apply(ave, 1, scale_fc))


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
dr = as.matrix(ave_fc)
fit = kmeans(dr, 10)


dr_kmeans = log2(rep[order(fit$cluster),]+1)
dr_kmeans_plot = rbind(dr_kmeans, rep_nopk[used_id_nopk,])
dr_kmeans_plot[dr_kmeans_plot>7]=7
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.TS.', toString(nr), '.pdf', sep=''))
pheatmap(dr_kmeans_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

dr_kmeans_ave = log2(ave[order(fit$cluster),]+1)
dr_kmeans_ave_plot = rbind(dr_kmeans_ave, ave_nopk[used_id_nopk,])
dr_kmeans_ave_plot[dr_kmeans_ave_plot>7]=7
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('kmean.TS.ave.', toString(nr), '.pdf', sep=''))
pheatmap(dr_kmeans_ave_plot, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()








