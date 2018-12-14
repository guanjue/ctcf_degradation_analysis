### get parameters
args = commandArgs(trailingOnly=TRUE)

sample_list_file = args[1]
sample_list_file = 'id_list.txt'
sample_list = read.table(sample_list_file, header=FALSE)

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

d_chipseq_1_mat = c()
d_chipseq_2_mat = c()
d_chipseq_input_mat = c()

d_nonpk_1_mat = c()
d_nonpk_2_mat = c()
d_nonpk_input_mat = c()

output_vec = c()

### get signal matrix
for (i in c(1:4)){
	rep1_i = 1 + 2*(i-1)
	rep2_i = 2 + 2*(i-1)
	print(sample_list[rep1_i, 3])
	rep1 = paste(sample_list[rep1_i, 1], '.qtPCR_regions.sorted.rc_per_bp.tab', sep='')
	rep2 = paste(sample_list[rep2_i, 1], '.qtPCR_regions.sorted.rc_per_bp.tab', sep='')
	input = paste(sample_list[rep1_i, 2], '.qtPCR_regions.input_sig.macs.tab', sep='')
	rep1_nonpk = paste(sample_list[rep1_i, 1], '.mm9_500_nonpk_r.sorted.rc_per_bp.tab', sep='')
	rep2_nonpk = paste(sample_list[rep2_i, 1], '.mm9_500_nonpk_r.sorted.rc_per_bp.tab', sep='')
	input_nonpk = paste(sample_list[rep1_i, 2], '.mm9_500_nonpk_r.input_sig.macs.tab', sep='')
	qtpcr = paste('qtPCR_sig.', sample_list[rep1_i, 3], '.txt', sep='')
	output = paste(sample_list[rep1_i, 3], '.chipseq_vs_qtpcr.png', sep='')
	output_vec[i] = output
	### read input files
	d_qtpcr = scan(qtpcr)*1000
	d_chipseq_1 = scan(rep1)*1000
	d_chipseq_2 = scan(rep2)*1000
	d_chipseq_input = scan(input)*1000
	d_chipseq_1_nonpk = scan(rep1_nonpk)*1000
	d_chipseq_2_nonpk = scan(rep2_nonpk)*1000
	d_chipseq_input_nonpk = scan(input_nonpk)*1000
	### get signal matrix
	d_qtpcr_mat = cbind(d_qtpcr_mat, d_qtpcr)
	d_chipseq_1_mat = cbind(d_chipseq_1_mat, d_chipseq_1)
	d_chipseq_2_mat = cbind(d_chipseq_2_mat, d_chipseq_2)
	d_chipseq_input_mat = cbind(d_chipseq_input_mat, d_chipseq_input)
	d_nonpk_1_mat = cbind(d_nonpk_1_mat, d_chipseq_1_nonpk)
	d_nonpk_2_mat = cbind(d_nonpk_2_mat, d_chipseq_2_nonpk)
	d_nonpk_input_mat = cbind(d_nonpk_input_mat, d_chipseq_input_nonpk)
}

### plot ChIP-seq replicates signal at qtPCR regions
lim_min_chipseq = min(cbind(d_chipseq_1_mat, d_chipseq_2_mat))
lim_max_chipseq = max(cbind(d_chipseq_1_mat, d_chipseq_2_mat))
color_list = c('black', 'purple', 'blue', 'green', 'orange', 'red')
pdf('overall_chipseq_rep_scatter_plot.pdf')
plot(d_chipseq_1_mat[,1], d_chipseq_2_mat[,1], pch=16, col=color_list[1], xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq))
for (i in c(2:6)){
	points(d_chipseq_1_mat[,i], d_chipseq_2_mat[,i], pch=16, col=color_list[i])
}
abline(0, 1, lwd=1.5, col = 'red')
dev.off()

### extract overall scale model for qtPCR
#( d_qtpcr_mat - mean(d_qtpcr_mat) ) / sd(d_qtpcr_mat) ~ ( d_chipseq_i_mat - mean(d_chipseq_i_mat) ) / sd(d_chipseq_i_mat)
d_chipseq_1_mat_fc = d_chipseq_1_mat / d_chipseq_input_mat
d_chipseq_2_mat_fc = d_chipseq_2_mat / d_chipseq_input_mat

d_qtpcr_mat_scale = ( d_qtpcr_mat - mean(d_qtpcr_mat) ) / sd(d_qtpcr_mat) * sd(cbind(d_chipseq_1_mat_fc, d_chipseq_2_mat_fc)) + mean(cbind(d_chipseq_1_mat_fc, d_chipseq_2_mat_fc)) 

### plot qtPCR replicates signal at qtPCR regions
lim_min_chipseq = min(cbind((d_chipseq_1_mat_fc), (d_chipseq_2_mat_fc)))
lim_max_chipseq = max(cbind((d_chipseq_1_mat_fc), (d_chipseq_2_mat_fc)))

pdf('overall_chipseq_vs_qtpcr_scatter_plot.pdf')
plot(d_qtpcr_mat_scale, d_chipseq_1_mat/d_chipseq_input_mat, pch=16, col='blue', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), log='')
points(d_qtpcr_mat_scale, d_chipseq_2_mat/d_chipseq_input_mat, pch=16, col='red')
#points(d_qtpcr_mat_scale, (d_chipseq_1_mat[,i]+d_chipseq_2_mat)/2/d_chipseq_input_mat, pch=16, col='black')
abline(0, 1, lwd=1.5, col = 'red')
dev.off()

pdf('overall_chipseq_vs_qtpcr_scatter_plot_peaktype_col.pdf')
plot(cbind(d_qtpcr_mat_scale,d_qtpcr_mat_scale)[tolower(peak_type)=='strong-binding',], cbind(d_chipseq_1_mat/d_chipseq_input_mat,d_chipseq_2_mat/d_chipseq_input_mat)[tolower(peak_type)=='strong-binding',], pch=16, col='red', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), log='')
points(cbind(d_qtpcr_mat_scale,d_qtpcr_mat_scale)[tolower(peak_type)=='medium-binding',], cbind(d_chipseq_1_mat/d_chipseq_input_mat,d_chipseq_2_mat/d_chipseq_input_mat)[tolower(peak_type)=='medium-binding',], pch=16, col='green')
points(cbind(d_qtpcr_mat_scale,d_qtpcr_mat_scale)[tolower(peak_type)=='weak-binding',], cbind(d_chipseq_1_mat/d_chipseq_input_mat,d_chipseq_2_mat/d_chipseq_input_mat)[tolower(peak_type)=='weak-binding',], pch=16, col='blue')
points(mean(cbind(d_qtpcr_mat_scale,d_qtpcr_mat_scale)[tolower(peak_type)=='strong-binding',]), mean(cbind(d_chipseq_1_mat/d_chipseq_input_mat,d_chipseq_2_mat/d_chipseq_input_mat)[tolower(peak_type)=='strong-binding',]), pch=15, cex=1.5, col='red')
points(mean(cbind(d_qtpcr_mat_scale,d_qtpcr_mat_scale)[tolower(peak_type)=='medium-binding',]), mean(cbind(d_chipseq_1_mat/d_chipseq_input_mat,d_chipseq_2_mat/d_chipseq_input_mat)[tolower(peak_type)=='medium-binding',]), pch=15, cex=1.5, col='green')
points(mean(cbind(d_qtpcr_mat_scale,d_qtpcr_mat_scale)[tolower(peak_type)=='weak-binding',]), mean(cbind(d_chipseq_1_mat/d_chipseq_input_mat,d_chipseq_2_mat/d_chipseq_input_mat)[tolower(peak_type)=='weak-binding',]), pch=15, cex=1.5, col='blue')
#points(d_qtpcr_mat_scale, (d_chipseq_1_mat[,i]+d_chipseq_2_mat)/2/d_chipseq_input_mat, pch=16, col='black')
abline(0, 1, lwd=1.5, col = 'red')
dev.off()


AB_1_mat = c()
AB_2_mat = c()
AB_average_mat = c()

AB_TS1_mat = c()
AB_TS2_mat = c()
AB_TSaverage_mat = c()


### S3norm vs TSnorm and plot
for (i in c(1:4)){
	print(i)
	smallnum = 1
	output_i = output_vec[i]
	d_qtpcr_mat_scale_i = d_qtpcr_mat_scale[,i]
	d_chipseq_1_mat_FC_i = (d_chipseq_1_mat[,i]+smallnum)/(d_chipseq_input_mat[,i]+smallnum)
	d_chipseq_2_mat_FC_i = (d_chipseq_2_mat[,i]+smallnum)/(d_chipseq_input_mat[,i]+smallnum)
	d_chipseq_average_mat_FC_i = ((d_chipseq_1_mat[,i]+d_chipseq_2_mat[,i])/2+smallnum)/(d_chipseq_input_mat[,i]+smallnum)
	d_chipseq_1_mat_FC_i_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_1_mat_FC_i)
	d_chipseq_2_mat_FC_i_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_2_mat_FC_i)
	d_chipseq_average_mat_FC_i_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i)

	d_nonpk_1_mat_FC_i_tar = (d_nonpk_1_mat[,i]+smallnum)/(d_nonpk_input_mat[,i]+smallnum)
	d_nonpk_2_mat_FC_i_tar = (d_nonpk_2_mat[,i]+smallnum)/(d_nonpk_input_mat[,i]+smallnum)
	d_nonpk_average_mat_FC_i_tar = ((d_nonpk_1_mat[,i]+d_nonpk_2_mat[,i])/2+smallnum)/(d_nonpk_input_mat[,i]+smallnum)
	d_nonpk_1_mat_FC_i_ref = (d_nonpk_1_mat[,4]+smallnum)/(d_nonpk_input_mat[,4]+smallnum)
	d_nonpk_2_mat_FC_i_ref = (d_nonpk_2_mat[,4]+smallnum)/(d_nonpk_input_mat[,4]+smallnum)
	d_nonpk_average_mat_FC_i_ref = ((d_nonpk_1_mat[,4]+d_nonpk_2_mat[,4])/2+smallnum)/(d_nonpk_input_mat[,4]+smallnum)

	#pdf(output_i, width=14, height=14)
	png(output_i, width=800, height=800)
	print(summary(d_nonpk_1_mat_FC_i_tar))
	par(mfrow=c(2,2))
	plot(d_nonpk_1_mat_FC_i_ref, d_nonpk_1_mat_FC_i_tar, pch=16, col='gray', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), main = paste('R2_rep1:', round(d_chipseq_1_mat_FC_i_cor, 3), 'R2_rep2:', round(d_chipseq_2_mat_FC_i_cor, 3), 'R2_average:', round(d_chipseq_average_mat_FC_i_cor, 3)))
	points(d_qtpcr_mat_scale_i, d_chipseq_1_mat_FC_i, pch=16, col='purple')
	points(d_qtpcr_mat_scale_i, d_chipseq_2_mat_FC_i, pch=16, col='brown')
	points(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i, pch=16, col='black')
	points(mean(d_qtpcr_mat_scale_i), mean(d_chipseq_average_mat_FC_i), pch=16, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='strong-binding']), mean(d_chipseq_1_mat_FC_i[tolower(peak_type)=='strong-binding']), pch=15, cex = 1.5, col='red')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='medium-binding']), mean(d_chipseq_2_mat_FC_i[tolower(peak_type)=='medium-binding']), pch=15, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='weak-binding']), mean(d_chipseq_average_mat_FC_i[tolower(peak_type)=='weak-binding']), pch=15, cex = 1.5, col='blue')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)!='strong-binding']), mean(d_chipseq_average_mat_FC_i[tolower(peak_type)!='strong-binding']), pch=15, cex = 1.5, col='black')
	abline(0, 1, lwd=1.5, col = 'red')
	TSnorm1 = 1.0 / mean(d_chipseq_1_mat_FC_i) * mean(d_qtpcr_mat_scale_i)
	TSnorm2 = 1.0 / mean(d_chipseq_2_mat_FC_i) * mean(d_qtpcr_mat_scale_i)
	TSnorm_average = 1.0 / mean(d_chipseq_average_mat_FC_i) * mean(d_qtpcr_mat_scale_i)
	d_chipseq_1_mat_FC_i_scale = d_chipseq_1_mat_FC_i * TSnorm1
	d_chipseq_2_mat_FC_i_scale = d_chipseq_2_mat_FC_i * TSnorm2
	d_nonpk_1_mat_FC_i_tar_scale = d_nonpk_1_mat_FC_i_tar * TSnorm1
	d_chipseq_average_mat_FC_i_scale = d_chipseq_average_mat_FC_i * TSnorm_average
	d_chipseq_1_mat_FC_i_scale_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_1_mat_FC_i_scale)
	d_chipseq_2_mat_FC_i_scale_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_2_mat_FC_i_scale)
	d_chipseq_average_mat_FC_i_scale_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i_scale)
	plot(d_nonpk_1_mat_FC_i_ref, d_nonpk_1_mat_FC_i_tar_scale, pch=16, col='gray', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), main = paste('R2_rep1:', round(d_chipseq_1_mat_FC_i_cor, 3), 'R2_rep2:', round(d_chipseq_2_mat_FC_i_cor, 3), 'R2_average:', round(d_chipseq_average_mat_FC_i_cor, 3)))
	points(d_qtpcr_mat_scale_i, d_chipseq_1_mat_FC_i_scale, pch=16, col='purple')
	points(d_qtpcr_mat_scale_i, d_chipseq_2_mat_FC_i_scale, pch=16, col='brown')
	points(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i_scale, pch=16, col='black')
	points(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i_scale, pch=16, col='black')
	points(mean(d_qtpcr_mat_scale_i), mean(d_chipseq_average_mat_FC_i_scale), pch=16, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='strong-binding']), mean(d_chipseq_average_mat_FC_i_scale[tolower(peak_type)=='strong-binding']), pch=15, cex = 1.5, col='red')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='medium-binding']), mean(d_chipseq_average_mat_FC_i_scale[tolower(peak_type)=='medium-binding']), pch=15, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='weak-binding']), mean(d_chipseq_average_mat_FC_i_scale[tolower(peak_type)=='weak-binding']), pch=15, cex = 1.5, col='blue')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)!='strong-binding']), mean(d_chipseq_average_mat_FC_i_scale[tolower(peak_type)!='strong-binding']), pch=15, cex = 1.5, col='black')
	abline(0, 1, lwd=1.5, col = 'red')
	### s3norm
	d_chipseq_1_mat_FC_i_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_1_mat_FC_i)
	d_chipseq_2_mat_FC_i_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_2_mat_FC_i)
	d_chipseq_average_mat_FC_i_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i)
	plot(d_qtpcr_mat_scale_i, d_chipseq_1_mat_FC_i, pch=16, col='purple', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), main = paste('R2_rep1:', round(d_chipseq_1_mat_FC_i_cor, 3), 'R2_rep2:', round(d_chipseq_2_mat_FC_i_cor, 3), 'R2_average:', round(d_chipseq_average_mat_FC_i_cor, 3)))
	points(d_qtpcr_mat_scale_i, d_chipseq_2_mat_FC_i, pch=16, col='brown')
	points(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i, pch=16, col='black')
	points(mean(d_qtpcr_mat_scale_i), mean(d_chipseq_average_mat_FC_i), pch=16, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='strong-binding']), mean(d_chipseq_1_mat_FC_i[tolower(peak_type)=='strong-binding']), pch=15, cex = 1.5, col='red')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='medium-binding']), mean(d_chipseq_2_mat_FC_i[tolower(peak_type)=='medium-binding']), pch=15, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='weak-binding']), mean(d_chipseq_average_mat_FC_i[tolower(peak_type)=='weak-binding']), pch=15, cex = 1.5, col='blue')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)!='strong-binding']), mean(d_chipseq_average_mat_FC_i[tolower(peak_type)!='strong-binding']), pch=15, cex = 1.5, col='black')
	abline(0, 1, lwd=1.5, col = 'red')
	used_pk = peak_type!='1Strong-Stable'
	#AB_1 = s3norm(d_qtpcr_mat_scale_i, d_chipseq_1_mat_FC_i, peak_type=='Strong-Stable', peak_type!='Strong-Stable', 1.0, 2.0, 1e-5)
	AB_1 = s3norm_nonpkbg(d_qtpcr_mat_scale_i[used_pk], d_chipseq_1_mat_FC_i[used_pk], d_nonpk_1_mat_FC_i_ref, d_nonpk_1_mat_FC_i_tar, 1.0, 2.0, 1e-5)
	#AB_1 = s3norm_nonpkbg_matchmean(log2(d_qtpcr_mat_scale_i), log2(d_chipseq_1_mat_FC_i), log2(d_nonpk_1_mat_FC_i_ref), log2(d_nonpk_1_mat_FC_i_tar))
	d_chipseq_1_mat_FC_i_s3norm = AB_1[1] * (d_chipseq_1_mat_FC_i ^ AB_1[2])
	d_nonpk_1_mat_FC_i_tar_s3norm = AB_1[1] * (d_nonpk_1_mat_FC_i_tar ^ AB_1[2])
	#AB_2 = s3norm(d_qtpcr_mat_scale_i, d_chipseq_2_mat_FC_i, peak_type=='Strong-Stable', peak_type!='Strong-Stable', 1.0, 2.0, 1e-5)
	AB_2 = s3norm_nonpkbg(d_qtpcr_mat_scale_i[used_pk], d_chipseq_2_mat_FC_i[used_pk], d_nonpk_2_mat_FC_i_ref, d_nonpk_2_mat_FC_i_tar, 1.0, 2.0, 1e-5)
	#AB_2 = s3norm_nonpkbg_matchmean(log2(d_qtpcr_mat_scale_i), log2(d_chipseq_2_mat_FC_i), log2(d_nonpk_2_mat_FC_i_ref), log2(d_nonpk_2_mat_FC_i_tar))
	d_chipseq_2_mat_FC_i_s3norm = AB_2[1] * (d_chipseq_2_mat_FC_i ^ AB_2[2])
	#AB_average = s3norm(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i, peak_type=='Strong-Stable', peak_type!='Strong-Stable', 1.0, 2.0, 1e-5)
	AB_average = s3norm_nonpkbg(d_qtpcr_mat_scale_i[used_pk], d_chipseq_average_mat_FC_i[used_pk], d_nonpk_average_mat_FC_i_ref, d_nonpk_average_mat_FC_i_tar, 1.0, 2.0, 1e-5)
	#AB_average = s3norm_nonpkbg_matchmean(log2(d_qtpcr_mat_scale_i), log2(d_chipseq_average_mat_FC_i), log2(d_nonpk_average_mat_FC_i_ref), log2(d_nonpk_average_mat_FC_i_tar))
	d_chipseq_average_mat_FC_i_s3norm = AB_average[1] * (d_chipseq_average_mat_FC_i ^ AB_average[2])
	d_chipseq_1_mat_FC_i_s3norm_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_1_mat_FC_i_s3norm)
	d_chipseq_2_mat_FC_i_s3norm_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_2_mat_FC_i_s3norm)
	d_chipseq_average_mat_FC_i_s3norm_cor = cor(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i_s3norm)
	plot(d_nonpk_1_mat_FC_i_ref, d_nonpk_1_mat_FC_i_tar_s3norm, pch=16, col='gray', xlim=c(lim_min_chipseq, lim_max_chipseq), ylim=c(lim_min_chipseq, lim_max_chipseq), main = paste('R2_rep1:', round(d_chipseq_1_mat_FC_i_s3norm_cor, 3), 'R2_rep2:', round(d_chipseq_2_mat_FC_i_s3norm_cor, 3), 'R2_average:', round(d_chipseq_average_mat_FC_i_s3norm_cor, 3)))
	points(d_qtpcr_mat_scale_i, d_chipseq_1_mat_FC_i_s3norm, pch=16, col='purple')
	points(d_qtpcr_mat_scale_i, d_chipseq_2_mat_FC_i_s3norm, pch=16, col='brown')
	points(d_qtpcr_mat_scale_i, d_chipseq_average_mat_FC_i_s3norm, pch=16, col='black')
	points(mean(d_qtpcr_mat_scale_i), mean(d_chipseq_average_mat_FC_i_s3norm), pch=16, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='strong-binding']), mean(d_chipseq_1_mat_FC_i_s3norm[tolower(peak_type)=='strong-binding']), pch=15, cex = 1.5, col='red')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='medium-binding']), mean(d_chipseq_2_mat_FC_i_s3norm[tolower(peak_type)=='medium-binding']), pch=15, cex = 1.5, col='green')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='weak-binding']), mean(d_chipseq_average_mat_FC_i_s3norm[tolower(peak_type)=='weak-binding']), pch=15, cex = 1.5, col='blue')
	points(mean(d_qtpcr_mat_scale_i[tolower(peak_type)!='strong-binding']), mean(d_chipseq_average_mat_FC_i_s3norm[tolower(peak_type)!='strong-binding']), pch=15, cex = 1.5, col='black')
	abline(0, 1, lwd=1.5, col = 'red')
	dev.off()
	AB_1_mat = rbind(AB_1_mat, AB_1)
	AB_2_mat = rbind(AB_2_mat, AB_2)
	AB_average_mat = rbind(AB_average_mat, AB_average)
	AB_TS1_mat = rbind(AB_TS1_mat, c(TSnorm1, 1))
	AB_TS2_mat = rbind(AB_TS2_mat, c(TSnorm2, 1))
	AB_TSaverage_mat = rbind(AB_TSaverage_mat, c(TSnorm_average, 1))
	print(mean(d_qtpcr_mat_scale_i[tolower(peak_type)=='strong-binding']))
	print(mean(d_qtpcr_mat_scale_i[tolower(peak_type)!='strong-binding']))
	print(mean(d_chipseq_average_mat_FC_i_s3norm[tolower(peak_type)=='strong-binding']))
	print(mean(d_chipseq_average_mat_FC_i_s3norm[tolower(peak_type)!='strong-binding']))

}






	print(mean(d_qtpcr_mat_scale_i[peak_type=='Strong-Stable']))
	print(mean(d_qtpcr_mat_scale_i[peak_type!='Strong-Stable']))
	print(mean(d_chipseq_1_mat_FC_i_s3norm[peak_type=='Strong-Stable']))
	print(mean(d_chipseq_1_mat_FC_i_s3norm[peak_type!='Strong-Stable']))
	print(mean(d_chipseq_2_mat_FC_i_s3norm[peak_type=='Strong-Stable']))
	print(mean(d_chipseq_2_mat_FC_i_s3norm[peak_type!='Strong-Stable']))

colnames(AB_1_mat) = c('A', 'B')
rownames(AB_1_mat) = sample_list[seq(1,8,2),3]
write.table(AB_1_mat, 'rep1_s3norm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

colnames(AB_2_mat) = c('A', 'B')
rownames(AB_2_mat) = sample_list[seq(1,8,2),3]
write.table(AB_2_mat, 'rep2_s3norm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

colnames(AB_average_mat) = c('A', 'B')
rownames(AB_average_mat) = sample_list[seq(1,8,2),3]
write.table(AB_average_mat, 'average_s3norm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

colnames(AB_TS1_mat) = c('A', 'B')
rownames(AB_TS1_mat) = sample_list[seq(1,8,2),3]
write.table(AB_TS1_mat, 'rep1_TSnorm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

colnames(AB_TS2_mat) = c('A', 'B')
rownames(AB_TS2_mat) = sample_list[seq(1,8,2),3]
write.table(AB_TS2_mat, 'rep2_TSnorm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)

colnames(AB_TSaverage_mat) = c('A', 'B')
rownames(AB_TSaverage_mat) = sample_list[seq(1,8,2),3]
write.table(AB_TSaverage_mat, 'average_TSnorm_AB.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)


