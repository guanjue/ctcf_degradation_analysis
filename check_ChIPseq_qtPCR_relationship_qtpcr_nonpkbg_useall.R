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
		d_chipseq_1 = scan(rep1)#*1000
		d_chipseq_mat = cbind(d_chipseq_mat, d_chipseq_1)
		d_chipseq_1_nonpk = scan(rep1_nonpk)#*1000
		d_nonpk_mat = cbind(d_nonpk_mat, d_chipseq_1_nonpk)
		input = paste(input_list, '.qtPCR_regions.input_sig.macs.tab', sep='')
		input_nonpk = paste(input_list, '.mm9_500.input_sig.macs.tab', sep='')
		qtpcr = paste('qtPCR_sig.', tp_type[i], '.txt', sep='')
		output = paste(tp_type[i], '.chipseq_vs_qtpcr.png', sep='')
		output_vec[i] = output
		### read input files
		d_qtpcr = scan(qtpcr)#*1000
		d_chipseq_input = scan(input)#*1000
		d_chipseq_input_nonpk = scan(input_nonpk)#*1000
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
d_chipseq_mat_fc = d_chipseq_mat #/ d_chipseq_input_mat

d_qtpcr_mat_scale = ( d_qtpcr_mat - mean(d_qtpcr_mat) ) / sd(d_qtpcr_mat) * sd(d_chipseq_mat_fc) + mean(d_chipseq_mat_fc) 

### plot qtPCR replicates signal at qtPCR regions
lim_min_chipseq = min(d_chipseq_mat_fc)
lim_max_chipseq = max(d_chipseq_mat_fc)





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

	d_qtpcr_mat_scale_i = d_qtpcr_mat_scale[,i]

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
	d_nonpk_average_mat_FC_i_ref = (rowMeans(d_nonpk_mat[,used_colid_ref])+smallnum)/(rowMeans(d_nonpk_input_mat[,used_colid_ref])+smallnum)

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


