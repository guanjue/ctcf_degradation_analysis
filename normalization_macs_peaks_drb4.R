### get parameters
args = commandArgs(trailingOnly=TRUE)

sample_list_file = args[1]
sample_list_file = 'id_list.txt'
sample_list = read.table(sample_list_file, header=FALSE)

### get pk num and get sample id
file1 = paste(sample_list[1, 1], '.allmerged.sorted.rc_per_bp.tab', sep='')
file1 = scan(file1)
pknum = length(file1)
set.seed(2018)
used_id = sample(pknum, 10000)

peak_type = read.table('qtPCR_regions.sorted.bed.txt', header=FALSE)[,5]

### read S3norm AB
AB_1_used = read.table('rep1_TSnorm_AB.txt', sep='\t', header=TRUE)
AB_2_used = read.table('rep2_TSnorm_AB.txt', sep='\t', header=TRUE)
AB_average_used = read.table('average_TSnorm_AB.txt', sep='\t', header=TRUE)

AB_1_used = read.table('rep1_s3norm_AB.txt', sep='\t', header=TRUE)
AB_2_used = read.table('rep2_s3norm_AB.txt', sep='\t', header=TRUE)
AB_average_used = read.table('average_s3norm_AB.txt', sep='\t', header=TRUE)

rep1_sig_s3norm_mat = c()
rep2_sig_s3norm_mat = c()
rep1_sig_TSnorm_mat = c()
rep2_sig_TSnorm_mat = c()
rep1_sig_s3norm_mat_nopk = c()
rep2_sig_s3norm_mat_nopk = c()

### normalization
for (i in c(1:4)){
	rep1_i = 1 + 2*(i-1)
	rep2_i = 2 + 2*(i-1)
	output = sample_list[rep1_i, 3]
	print(output)
	### read macs peaks
	rep1 = paste(sample_list[rep1_i, 1], '.allmerged.sorted.rc_per_bp.tab', sep='')
	rep2 = paste(sample_list[rep2_i, 1], '.allmerged.sorted.rc_per_bp.tab', sep='')
	input = paste(sample_list[rep1_i, 2], '.allmerged.input_sig.macs.tab', sep='')
	### read data
	rep1_sig = scan(rep1)*1000
	rep2_sig = scan(rep2)*1000
	input_sig = scan(input)*1000
	### read qtpcr
	rep1_qtpcr = paste(sample_list[rep1_i, 1], '.qtPCR_regions.sorted.rc_per_bp.tab', sep='')
	rep2_qtpcr = paste(sample_list[rep2_i, 1], '.qtPCR_regions.sorted.rc_per_bp.tab', sep='')
	input_qtpcr = paste(sample_list[rep1_i, 2], '.qtPCR_regions.input_sig.macs.tab', sep='')
	### read data
	rep1_sig_qtpcr = scan(rep1_qtpcr)*1000
	rep2_sig_qtpcr = scan(rep2_qtpcr)*1000
	input_sig_qtpcr = scan(input_qtpcr)*1000
	### read nopeak
	rep1_nopk = paste(sample_list[rep1_i, 1], '.mm9_500_nonpk_r.sorted.rc_per_bp.tab', sep='')
	rep2_nopk = paste(sample_list[rep2_i, 1], '.mm9_500_nonpk_r.sorted.rc_per_bp.tab', sep='')
	input_nopk = paste(sample_list[rep1_i, 2], '.mm9_500_nonpk_r.input_sig.macs.tab', sep='')
	### read data
	rep1_sig_nopk = scan(rep1_nopk)*1000
	rep2_sig_nopk = scan(rep2_nopk)*1000
	input_sig_nopk = scan(input_nopk)*1000

	### get S3norm A & B
	AB_1_i = AB_1_used[i,]
	AB_2_i = AB_2_used[i,]
	AB_average_i = AB_average_used[i,]
	### normalize
	smallnum = 1
	rep1_sig_fc = (rep1_sig+smallnum)/(input_sig+smallnum)
	rep2_sig_fc = (rep2_sig+smallnum)/(input_sig+smallnum)
	rep1_sig_s3norm = (rep1_sig_fc)^AB_1_i[1,2] * AB_1_i[1,1]
	rep2_sig_s3norm = (rep2_sig_fc)^AB_2_i[1,2] * AB_2_i[1,1]
	rep1_sig_s3norm_mat = cbind(rep1_sig_s3norm_mat, rep1_sig_s3norm)
	rep2_sig_s3norm_mat = cbind(rep2_sig_s3norm_mat, rep2_sig_s3norm)

	rep1_sig_TSnorm_mat = cbind(rep1_sig_TSnorm_mat, rep1_sig_fc/sum(rep1_sig_fc)*1000000)
	rep2_sig_TSnorm_mat = cbind(rep2_sig_TSnorm_mat, rep2_sig_fc/sum(rep2_sig_fc)*1000000)
	### qtpcr
	rep1_sig_fc_qtpcr = rep1_sig_qtpcr/input_sig_qtpcr
	rep2_sig_fc_qtpcr = rep2_sig_qtpcr/input_sig_qtpcr
	rep1_sig_s3norm_qtpcr = (rep1_sig_fc_qtpcr)^AB_1_i[1,2] * AB_1_i[1,1]
	rep2_sig_s3norm_qtpcr = (rep2_sig_fc_qtpcr)^AB_2_i[1,2] * AB_2_i[1,1]
	### nopk
	rep1_sig_fc_nopk = rep1_sig_nopk/input_sig_nopk
	rep2_sig_fc_nopk = rep2_sig_nopk/input_sig_nopk
	rep1_sig_s3norm_nopk = (rep1_sig_fc_nopk)^AB_1_i[1,2] * AB_1_i[1,1]
	rep2_sig_s3norm_nopk = (rep2_sig_fc_nopk)^AB_2_i[1,2] * AB_2_i[1,1]
	rep1_sig_s3norm_mat_nopk = cbind(rep1_sig_s3norm_mat_nopk, rep1_sig_s3norm_nopk)
	rep2_sig_s3norm_mat_nopk = cbind(rep2_sig_s3norm_mat_nopk, rep2_sig_s3norm_nopk)

	### plot replicates signal
	pdf(paste(output, '_s3norm_replicates.pdf', sep=''), width=7, height=14)
	par(mfrow=c(2,1))
	plot(rep1_sig_fc[used_id], rep2_sig_fc[used_id], xlim=c(0,150), ylim=c(0,150), main=round(cor(rep1_sig_fc[used_id], rep2_sig_fc[used_id]), 3))
	points((rep1_sig_fc_qtpcr[peak_type=='Strong-Stable']), (rep2_sig_fc_qtpcr[peak_type=='Strong-Stable']), pch=16, col='red')
	points((rep1_sig_fc_qtpcr[peak_type=='Strong-Dynamic']), (rep2_sig_fc_qtpcr[peak_type=='Strong-Dynamic']), pch=16, col='green')
	points((rep1_sig_fc_qtpcr[peak_type=='Weak']), (rep2_sig_fc_qtpcr[peak_type=='Weak']), pch=16, col='blue')
	points(mean(rep1_sig_fc_qtpcr[peak_type=='Strong-Stable']), mean(rep2_sig_fc_qtpcr[peak_type=='Strong-Stable']), pch=15, cex = 1.5, col='red')
	points(mean(rep1_sig_fc_qtpcr[peak_type=='Strong-Dynamic']), mean(rep2_sig_fc_qtpcr[peak_type=='Strong-Dynamic']), pch=15, cex = 1.5, col='green')
	points(mean(rep1_sig_fc_qtpcr[peak_type=='Weak']), mean(rep2_sig_fc_qtpcr[peak_type=='Weak']), pch=15, cex = 1.5, col='blue')
	points(mean(rep1_sig_fc_qtpcr[peak_type!='Strong-Stable']), mean(rep2_sig_fc_qtpcr[peak_type!='Strong-Stable']), pch=15, cex = 1.5, col='black')
	abline(0,1,col='red', lwd=1.5)
	plot(rep1_sig_s3norm[used_id], rep2_sig_s3norm[used_id], xlim=c(0,150), ylim=c(0,150), main=round(cor(rep1_sig_s3norm[used_id], rep2_sig_s3norm[used_id]), 3))
	points((rep1_sig_s3norm_qtpcr[peak_type=='Strong-Stable']), (rep2_sig_s3norm_qtpcr[peak_type=='Strong-Stable']), pch=16, col='red')
	points((rep1_sig_s3norm_qtpcr[peak_type=='Strong-Dynamic']), (rep2_sig_s3norm_qtpcr[peak_type=='Strong-Dynamic']), pch=16, col='green')
	points((rep1_sig_s3norm_qtpcr[peak_type=='Weak']), (rep2_sig_s3norm_qtpcr[peak_type=='Weak']), pch=16, col='blue')
	points(mean(rep1_sig_s3norm_qtpcr[peak_type=='Strong-Stable']), mean(rep2_sig_s3norm_qtpcr[peak_type=='Strong-Stable']), pch=15, cex = 1.5, col='red')
	points(mean(rep1_sig_s3norm_qtpcr[peak_type=='Strong-Dynamic']), mean(rep2_sig_s3norm_qtpcr[peak_type=='Strong-Dynamic']), pch=15, cex = 1.5, col='green')
	points(mean(rep1_sig_s3norm_qtpcr[peak_type=='Weak']), mean(rep2_sig_s3norm_qtpcr[peak_type=='Weak']), pch=15, cex = 1.5, col='blue')
	points(mean(rep1_sig_s3norm_qtpcr[peak_type!='Strong-Stable']), mean(rep2_sig_s3norm_qtpcr[peak_type!='Strong-Stable']), pch=15, cex = 1.5, col='black')
	abline(0,1,col='red', lwd=1.5)
	dev.off()
	print(mean(rep1_sig_s3norm_qtpcr[peak_type=='Strong-Stable']))
	print(mean(rep1_sig_s3norm_qtpcr[peak_type!='Strong-Stable']))
	print(mean(rep2_sig_s3norm_qtpcr[peak_type=='Strong-Stable']))
	print(mean(rep2_sig_s3norm_qtpcr[peak_type!='Strong-Stable']))
}


colnames(rep1_sig_s3norm_mat) = sample_list[seq(1,8,2),3]
colnames(rep2_sig_s3norm_mat) = sample_list[seq(1,8,2),3]

write.table(rep1_sig_s3norm_mat, 'rep1_sig_s3norm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep2_sig_s3norm_mat, 'rep2_sig_s3norm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(rep1_sig_s3norm_mat_nopk, 'rep1_sig_s3norm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep2_sig_s3norm_mat_nopk, 'rep2_sig_s3norm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table((rep1_sig_s3norm_mat+rep2_sig_s3norm_mat)/2, 'average_sig_s3norm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table((rep1_sig_s3norm_mat_nopk+rep2_sig_s3norm_mat_nopk)/2, 'average_sig_s3norm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)



write.table(rep1_sig_s3norm_mat, 'rep1_sig_TSnorm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep2_sig_s3norm_mat, 'rep2_sig_TSnorm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(rep1_sig_s3norm_mat_nopk, 'rep1_sig_TSnorm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep2_sig_s3norm_mat_nopk, 'rep2_sig_TSnorm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table((rep1_sig_s3norm_mat+rep2_sig_s3norm_mat)/2, 'average_sig_TSnorm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table((rep1_sig_s3norm_mat_nopk+rep2_sig_s3norm_mat_nopk)/2, 'average_sig_TSnorm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)












bed_file = read.table('allmerged.sorted.bed', header=FALSE)

rep1_sig_s3norm_mat_with_bed = cbind(bed_file[,4], as.data.frame(rep1_sig_s3norm_mat))[used_id,]
rep2_sig_s3norm_mat_with_bed = cbind(bed_file[,4], as.data.frame(rep2_sig_s3norm_mat))[used_id,]
colnames(rep1_sig_s3norm_mat_with_bed)[1] = 'bed'
colnames(rep2_sig_s3norm_mat_with_bed)[1] = 'bed'

write.table(rep1_sig_s3norm_mat_with_bed, 'rep1_sig_s3norm_mat.r.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep2_sig_s3norm_mat_with_bed, 'rep2_sig_s3norm_mat.r.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

smallnum = 1
rep1_sig_s3norm_mat_percentage = t(apply(rep1_sig_s3norm_mat, 1, function(x) (x+smallnum)/(x[1]+smallnum) ))
rep2_sig_s3norm_mat_percentage = t(apply(rep2_sig_s3norm_mat, 1, function(x) (x+smallnum)/(x[1]+smallnum) ))


rep1_sig_s3norm_mat_with_bed_percentage = cbind(bed_file[,4], as.data.frame(rep1_sig_s3norm_mat_percentage))[used_id,]
rep2_sig_s3norm_mat_with_bed_percentage = cbind(bed_file[,4], as.data.frame(rep2_sig_s3norm_mat_percentage))[used_id,]
colnames(rep1_sig_s3norm_mat_with_bed_percentage)[1] = 'bed'
colnames(rep2_sig_s3norm_mat_with_bed_percentage)[1] = 'bed'

write.table(rep1_sig_s3norm_mat_with_bed_percentage, 'rep1_sig_s3norm_mat_percentage.r.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep2_sig_s3norm_mat_with_bed_percentage, 'rep2_sig_s3norm_mat_percentage.r.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
--unscaled
--do_not_mean_center
s3norm_mat_merge = as.data.frame(bed_file[,4])
s3norm_mat_merge_colnames = as.data.frame(c('bed'))

for (i in c(1:6)){
	s3norm_mat_merge = cbind(s3norm_mat_merge, rep1_sig_s3norm_mat[,i])
	s3norm_mat_merge = cbind(s3norm_mat_merge, rep2_sig_s3norm_mat[,i])
	s3norm_mat_merge_colnames[2+(i-1)*2] = toString(sample_list[1+(i-1)*2,1])
	s3norm_mat_merge_colnames[3+(i-1)*2] = toString(sample_list[2+(i-1)*2,1])
}

colnames(s3norm_mat_merge) = s3norm_mat_merge_colnames
colnames(s3norm_mat_merge)[1] = 'bed'

write.table(s3norm_mat_merge, 's3norm_mat_merge.r.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)




### normalization
rep1_1 = 1 
rep2_1 = 2 
rep1_1 = paste(sample_list[1, 1], '.allmerged.sorted.rc_per_bp.tab', sep='')
rep2_1 = paste(sample_list[2, 1], '.allmerged.sorted.rc_per_bp.tab', sep='')
input_1 = paste(sample_list[1, 2], '.allmerged.input_sig.macs.tab', sep='')
rep1_qtpcr_1 = paste(sample_list[1, 1], '.qtPCR_regions.sorted.rc_per_bp.tab', sep='')
rep2_qtpcr_1 = paste(sample_list[2, 1], '.qtPCR_regions.sorted.rc_per_bp.tab', sep='')
input_qtpcr_1 = paste(sample_list[1, 2], '.qtPCR_regions.input_sig.macs.tab', sep='')

rep1_sig_1 = scan(rep1_1)
rep2_sig_1 = scan(rep2_1)
input_sig_1 = scan(input_1)
rep1_sig_qtpcr_1 = scan(rep1_qtpcr_1)
rep2_sig_qtpcr_1 = scan(rep2_qtpcr_1)
input_sig_qtpcr_1 = scan(input_qtpcr_1)
AB_1_i_1 = AB_1_used[1,]
AB_2_i_1 = AB_2_used[1,]
AB_average_i_1 = AB_average_used[1,]
rep1_sig_fc_1 = rep1_sig_1/input_sig_1
rep2_sig_fc_1 = rep2_sig_1/input_sig_1
rep1_sig_s3norm_1 = (rep1_sig_fc_1)^AB_1_i_1[1,2] * AB_1_i_1[1,1]
rep2_sig_s3norm_1 = (rep2_sig_fc_1)^AB_2_i_1[1,2] * AB_2_i_1[1,1]
### qtpcr
rep1_sig_fc_qtpcr_1 = rep1_sig_qtpcr_1/input_sig_qtpcr_1
rep2_sig_fc_qtpcr_1 = rep2_sig_qtpcr_1/input_sig_qtpcr_1
rep1_sig_s3norm_qtpcr_1 = (rep1_sig_fc_qtpcr_1)^AB_1_i_1[1,2] * AB_1_i_1[1,1]
rep2_sig_s3norm_qtpcr_1 = (rep2_sig_fc_qtpcr_1)^AB_2_i_1[1,2] * AB_2_i_1[1,1]
for (i in c(1:6)){
	rep1_i = 1 + 2*(i-1)
	rep2_i = 2 + 2*(i-1)
	output = sample_list[rep1_i, 3]
	print(output)
	### read macs peaks
	rep1 = paste(sample_list[rep1_i, 1], '.allmerged.sorted.rc_per_bp.tab', sep='')
	rep2 = paste(sample_list[rep2_i, 1], '.allmerged.sorted.rc_per_bp.tab', sep='')
	input = paste(sample_list[rep1_i, 2], '.allmerged.input_sig.macs.tab', sep='')
	### read data
	rep1_sig = scan(rep1)
	rep2_sig = scan(rep2)
	input_sig = scan(input)
	### read qtpcr
	rep1_qtpcr = paste(sample_list[rep1_i, 1], '.qtPCR_regions.sorted.rc_per_bp.tab', sep='')
	rep2_qtpcr = paste(sample_list[rep2_i, 1], '.qtPCR_regions.sorted.rc_per_bp.tab', sep='')
	input_qtpcr = paste(sample_list[rep1_i, 2], '.qtPCR_regions.input_sig.macs.tab', sep='')
	### read data
	rep1_sig_qtpcr = scan(rep1_qtpcr)
	rep2_sig_qtpcr = scan(rep2_qtpcr)
	input_sig_qtpcr = scan(input_qtpcr)
	### get S3norm A & B
	AB_1_i = AB_1_used[i,]
	AB_2_i = AB_2_used[i,]
	AB_average_i = AB_average_used[i,]
	### normalize
	rep1_sig_fc = rep1_sig/input_sig
	rep2_sig_fc = rep2_sig/input_sig
	rep1_sig_s3norm = (rep1_sig_fc)^AB_1_i[1,2] * AB_1_i[1,1]
	rep2_sig_s3norm = (rep2_sig_fc)^AB_2_i[1,2] * AB_2_i[1,1]
	### qtpcr
	rep1_sig_fc_qtpcr = rep1_sig_qtpcr/input_sig_qtpcr
	rep2_sig_fc_qtpcr = rep2_sig_qtpcr/input_sig_qtpcr
	rep1_sig_s3norm_qtpcr = (rep1_sig_fc_qtpcr)^AB_1_i[1,2] * AB_1_i[1,1]
	rep2_sig_s3norm_qtpcr = (rep2_sig_fc_qtpcr)^AB_2_i[1,2] * AB_2_i[1,1]
	### plot replicates signal
	#pdf(paste(output, '_norm_vs_diftime.pdf', sep=''), width=7, height=14)
	png(paste(output, '_s3norm_nonpk_vs_diftime.png', sep=''), width=400, height=800)
	par(mfrow=c(2,1))
	plot(rep1_sig_fc_1[used_id]/2+rep2_sig_fc_1[used_id]/2, rep1_sig_fc[used_id]/2+rep2_sig_fc[used_id]/2, xlim=c(0,150), ylim=c(0,150), main=round(cor(rep1_sig_fc[used_id], rep2_sig_fc[used_id]), 3))
	points((rep1_sig_fc_qtpcr_1[peak_type=='Strong-Stable'])/2+(rep2_sig_fc_qtpcr_1[peak_type=='Strong-Stable'])/2, (rep1_sig_fc_qtpcr[peak_type=='Strong-Stable'])/2+(rep2_sig_fc_qtpcr[peak_type=='Strong-Stable'])/2, pch=16, col='red')
	points((rep1_sig_fc_qtpcr_1[peak_type=='Strong-Dynamic'])/2+(rep2_sig_fc_qtpcr_1[peak_type=='Strong-Dynamic'])/2, (rep1_sig_fc_qtpcr[peak_type=='Strong-Dynamic'])/2+(rep2_sig_fc_qtpcr[peak_type=='Strong-Dynamic'])/2, pch=16, col='green')
	points((rep1_sig_fc_qtpcr_1[peak_type=='Weak'])/2+(rep2_sig_fc_qtpcr_1[peak_type=='Weak'])/2, (rep1_sig_fc_qtpcr[peak_type=='Weak'])/2+(rep2_sig_fc_qtpcr[peak_type=='Weak'])/2, pch=16, col='blue')
	points(mean(rep1_sig_fc_qtpcr_1[peak_type=='Strong-Stable'])/2+mean(rep2_sig_fc_qtpcr_1[peak_type=='Strong-Stable'])/2, mean(rep1_sig_fc_qtpcr[peak_type=='Strong-Stable'])/2+mean(rep2_sig_fc_qtpcr[peak_type=='Strong-Stable'])/2, pch=15, cex = 1.5, col='red')
	points(mean(rep1_sig_fc_qtpcr_1[peak_type=='Strong-Dynamic'])/2+mean(rep2_sig_fc_qtpcr_1[peak_type=='Strong-Dynamic'])/2, mean(rep1_sig_fc_qtpcr[peak_type=='Strong-Dynamic'])/2+mean(rep2_sig_fc_qtpcr[peak_type=='Strong-Dynamic'])/2, pch=15, cex = 1.5, col='green')
	points(mean(rep1_sig_fc_qtpcr_1[peak_type=='Weak'])/2+mean(rep2_sig_fc_qtpcr_1[peak_type=='Weak'])/2, mean(rep1_sig_fc_qtpcr[peak_type=='Weak'])/2+mean(rep2_sig_fc_qtpcr[peak_type=='Weak'])/2, pch=15, cex = 1.5, col='blue')
	points(mean(rep1_sig_fc_qtpcr_1[peak_type!='Strong-Stable'])/2+mean(rep2_sig_fc_qtpcr_1[peak_type!='Strong-Stable'])/2, mean(rep1_sig_fc_qtpcr[peak_type!='Strong-Stable'])/2+mean(rep2_sig_fc_qtpcr[peak_type!='Strong-Stable'])/2, pch=15, cex = 1.5, col='black')
	abline(0,1,col='red', lwd=1.5)
	plot(rep1_sig_s3norm_1[used_id]/2+rep2_sig_s3norm_1[used_id]/2, rep1_sig_s3norm[used_id]/2+rep2_sig_s3norm[used_id]/2, xlim=c(0,150), ylim=c(0,150), main=round(cor(rep1_sig_s3norm[used_id], rep2_sig_s3norm[used_id]), 3))
	points((rep1_sig_s3norm_qtpcr_1[peak_type=='Strong-Stable'])/2+(rep2_sig_s3norm_qtpcr_1[peak_type=='Strong-Stable'])/2, (rep1_sig_s3norm_qtpcr[peak_type=='Strong-Stable'])/2+(rep2_sig_s3norm_qtpcr[peak_type=='Strong-Stable'])/2, pch=16, col='red')
	points((rep1_sig_s3norm_qtpcr_1[peak_type=='Strong-Dynamic'])/2+(rep2_sig_s3norm_qtpcr_1[peak_type=='Strong-Dynamic'])/2, (rep1_sig_s3norm_qtpcr[peak_type=='Strong-Dynamic'])/2+(rep2_sig_s3norm_qtpcr[peak_type=='Strong-Dynamic'])/2, pch=16, col='green')
	points((rep1_sig_s3norm_qtpcr_1[peak_type=='Weak'])/2+(rep2_sig_s3norm_qtpcr_1[peak_type=='Weak'])/2, (rep1_sig_s3norm_qtpcr[peak_type=='Weak'])/2+(rep2_sig_s3norm_qtpcr[peak_type=='Weak'])/2, pch=16, col='blue')
	points(mean(rep1_sig_s3norm_qtpcr_1[peak_type=='Strong-Stable'])/2+mean(rep2_sig_s3norm_qtpcr_1[peak_type=='Strong-Stable'])/2, mean(rep1_sig_s3norm_qtpcr[peak_type=='Strong-Stable'])/2+mean(rep2_sig_s3norm_qtpcr[peak_type=='Strong-Stable'])/2, pch=15, cex = 1.5, col='red')
	points(mean(rep1_sig_s3norm_qtpcr_1[peak_type=='Strong-Dynamic'])/2+mean(rep2_sig_s3norm_qtpcr_1[peak_type=='Strong-Dynamic'])/2, mean(rep1_sig_s3norm_qtpcr[peak_type=='Strong-Dynamic'])/2+mean(rep2_sig_s3norm_qtpcr[peak_type=='Strong-Dynamic'])/2, pch=15, cex = 1.5, col='green')
	points(mean(rep1_sig_s3norm_qtpcr_1[peak_type=='Weak'])/2+mean(rep2_sig_s3norm_qtpcr_1[peak_type=='Weak'])/2, mean(rep1_sig_s3norm_qtpcr[peak_type=='Weak'])/2+mean(rep2_sig_s3norm_qtpcr[peak_type=='Weak'])/2, pch=15, cex = 1.5, col='blue')
	points(mean(rep1_sig_s3norm_qtpcr_1[peak_type!='Strong-Stable'])/2+mean(rep2_sig_s3norm_qtpcr_1[peak_type!='Strong-Stable'])/2, mean(rep1_sig_s3norm_qtpcr[peak_type!='Strong-Stable'])/2+mean(rep2_sig_s3norm_qtpcr[peak_type!='Strong-Stable'])/2, pch=15, cex = 1.5, col='black')
	abline(0,1,col='red', lwd=1.5)
	dev.off()
}





smallnum = 1
s3norm_mat_merge = read.table('rep1_sig_s3norm_mat.r.txt', header=T)
dr = as.matrix(s3norm_mat_merge[,-1])
rownames(dr) = s3norm_mat_merge[,1]

dr_fc = t(apply(dr, 1, function(x) (x+smallnum)/(mean(x)+smallnum)))
#dr_fc = dr
library(pheatmap)
library(mclust)

fit = kmeans(dr_fc, 7)
dr_kmeans = dr[order(fit$cluster),]
nr=10000
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
png(paste('kmean1.', toString(nr), '.png', sep=''))
pheatmap(dr_kmeans, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

BIC <- mclustBIC(dr_fc)
fit = Mclust(dr_fc, x = BIC)
dr_mclust = dr[order(fit$classification),]
nr=10000
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
png(paste('Mclust1.', toString(nr), '.png', sep=''))
pheatmap(dr_mclust, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()


s3norm_mat_merge = read.table('rep2_sig_s3norm_mat.r.txt', header=T)
dr = as.matrix(s3norm_mat_merge[,-1])
rownames(dr) = s3norm_mat_merge[,1]

dr_fc = t(apply(dr, 1, function(x) (x+smallnum)/(mean(x)+smallnum)))
#dr_fc = dr
library(pheatmap)
library(mclust)

fit = kmeans(dr_fc, 7)
dr_kmeans = dr[order(fit$cluster),]
nr=10000
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
png(paste('kmean2.', toString(nr), '.png', sep=''))
pheatmap(dr_kmeans, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

BIC <- mclustBIC(dr_fc)
fit = Mclust(dr_fc, x = BIC)
dr_mclust = dr[order(fit$classification),]
nr=10000
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
png(paste('Mclust2.', toString(nr), '.png', sep=''))
pheatmap(dr_mclust, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()





dr_fc = t(apply(dr, 1, function(x) (x+smallnum)/(mean(x)+smallnum)))
#dr_fc = dr
library(pheatmap)
library(mclust)

fit = kmeans(dr_fc, 7)
dr_kmeans = dr[order(fit$cluster),]
nr=10000
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
colnames(dr_kmeans) = c(1:12)
png(paste('kmean12.s3.', toString(nr), '.png', sep=''))
pheatmap(dr_kmeans, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=TRUE)
dev.off()

BIC <- mclustBIC(dr_fc)
fit = Mclust(dr_fc, x = BIC)
dr_mclust = dr[order(fit$classification),]
nr=10000
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
png(paste('Mclust12.', toString(nr), '.png', sep=''))
pheatmap(dr_mclust, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()



s3norm_mat_merge1 = read.table('rep1_sig_s3norm_mat.r.txt', header=T)
s3norm_mat_merge2 = read.table('rep2_sig_s3norm_mat.r.txt', header=T)

dr1 = as.matrix(s3norm_mat_merge1[,-1])
dr2 = as.matrix(s3norm_mat_merge2[,-1])

dr=c()
for (i in c(1:6)){
	dr = cbind(dr, dr1[,i])
	dr = cbind(dr, dr2[,i])
}

dr = (dr1+dr2)/2
rownames(dr) = s3norm_mat_merge1[,1]

dr_fc = t(apply(dr, 1, function(x) (x+smallnum)/(mean(x)+smallnum)))
#dr_fc = dr
library(pheatmap)
library(mclust)

fit = kmeans(dr_fc, 7)
dr_kmeans = dr[order(fit$cluster),]
nr=10000
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
colnames(dr_kmeans) = c(1:12)
png(paste('kmean12average.s3.', toString(nr), '.png', sep=''))
pheatmap(dr_kmeans, color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()


pdf('dif_0A_4A.pdf', width = 14, height = 14)
par(mfrow=c(2,2))
plot(rep1_sig_s3norm_mat[used_id, 1], rep1_sig_s3norm_mat[used_id, 2], xlim=c(0,150), ylim=c(0,150), main=round(cor(rep1_sig_fc[used_id], rep2_sig_fc[used_id]), 3))
abline(0, 1, col='red', lwd=1.5)
plot(rep2_sig_s3norm_mat[used_id, 1], rep2_sig_s3norm_mat[used_id, 2], xlim=c(0,150), ylim=c(0,150), main=round(cor(rep1_sig_fc[used_id], rep2_sig_fc[used_id]), 3))
abline(0, 1, col='red', lwd=1.5)
plot(rep1_sig_TSnorm_mat[used_id, 1], rep1_sig_TSnorm_mat[used_id, 2], xlim=c(0.1,150), ylim=c(0.1,150), main=round(cor(rep1_sig_fc[used_id], rep2_sig_fc[used_id]), 3), log='')
abline(0, 1, col='red', lwd=1.5)
plot(rep2_sig_TSnorm_mat[used_id, 1], rep2_sig_TSnorm_mat[used_id, 2], xlim=c(0.1,150), ylim=c(0.1,150), main=round(cor(rep1_sig_fc[used_id], rep2_sig_fc[used_id]), 3), log='')
abline(0, 1, col='red', lwd=1.5)
dev.off()











s3norm_mat_merge1 = read.table('rep1_sig_s3norm_mat.txt', header=T)
s3norm_mat_merge2 = read.table('rep2_sig_s3norm_mat.txt', header=T)

dr1 = as.matrix(s3norm_mat_merge1[,])
dr2 = as.matrix(s3norm_mat_merge2[,])

dr=c()
for (i in c(1:6)){
	dr = cbind(dr, dr1[,i])
	dr = cbind(dr, dr2[,i])
}

rownames(dr) = s3norm_mat_merge1[,1]

dr_average = (dr1+dr2)/2
rownames(dr_average) = s3norm_mat_merge1[,1]


smallnum=1
dr_average_fc = t(apply(dr_average, 1, function(x) (x+smallnum)/(mean(x)+smallnum)))

pdf('fc_hist.pdf', width = 7, height = 49)
par(mfrow=c(7,1))
hist(dr_average_fc, xlim=c(0,5), breaks=30)
abline(v=1, lwd=1.5, col='red')
box()
hist(dr_average_fc[,1], xlim=c(0,5), breaks=30)
abline(v=1, lwd=1.5, col='red')
box()
hist(dr_average_fc[,2], xlim=c(0,5), breaks=30)
abline(v=1, lwd=1.5, col='red')
box()
hist(dr_average_fc[,3], xlim=c(0,5), breaks=30)
abline(v=1, lwd=1.5, col='red')
box()
hist(dr_average_fc[,4], xlim=c(0,5), breaks=30)
abline(v=1, lwd=1.5, col='red')
box()
hist(dr_average_fc[,5], xlim=c(0,5), breaks=30)
abline(v=1, lwd=1.5, col='red')
box()
hist(dr_average_fc[,6], xlim=c(0,5), breaks=30)
abline(v=1, lwd=1.5, col='red')
box()
dev.off()


library(mclust)
dr_average_fc_vec = as.vector(dr_average_fc[,])
set.seed(2018)
mod_all_dr_average_fc_bic <- densityMclust(dr_average_fc_vec)
set.seed(2018)
mod_all_mod_all_dr_average_fc <- densityMclust(dr_average_fc_vec, G=3)

cluster_id = mod_all_mod_all_dr_average_fc$classification
cluster_mean = mod_all_mod_all_dr_average_fc$parameters$mean
print('2nd GMM cluster means: ')
print(cluster_mean)

rainbow_cp = rev(rainbow(length(cluster_mean)))


pdf('gmm.density.average_fc.pdf', width=14, height=7)
par(mfrow=c(1,2))
plot(mod_all_mod_all_dr_average_fc, what = "density", data = dr_average_fc_vec, breaks = 50)
#plotDensityMclust1(mod_all, data = signal_mat_log2_vec_high, hist.col = "lightgrey", hist.border = "white",  breaks = "Sturges", type = "persp")
for (i in c(1:length(cluster_mean))){
        print(i)
        x_input = seq(-4,9, 0.1)
        cp_i_mean = mod_all_mod_all_dr_average_fc$parameters$mean[i]
        cp_i_sd = (mod_all_mod_all_dr_average_fc$parameters$variance$sigmasq[i])^0.5
        cp_i_pro = mod_all_mod_all_dr_average_fc$parameters$pro[i]
        lines(x_input, cp_i_pro * dnorm(x_input, mean=cp_i_mean, sd=cp_i_sd), col=rainbow_cp[i])
}
plot(mod_all_dr_average_fc_bic, what = "BIC")
#plot(mod_all, what = "diagnostic", type = "cdf")
#plot(mod_all, what = "diagnostic", type = "qq")
dev.off()

signal_mat_index = dr_average_fc_vec
signal_thresh = min(signal_mat_index)

gmm_thresh = c()
gmm_thresh[1] = signal_thresh
for (i in c(2:(length(cluster_mean)))){
        print(i)
        gmm_thresh[i] = min(dr_average_fc_vec[cluster_id==i])
}
print(gmm_thresh)


dr_average_fc_binary = dr_average_fc
dr_average_fc_binary[dr_average_fc<=gmm_thresh[2]] = 0
dr_average_fc_binary[(dr_average_fc>gmm_thresh[2])*(dr_average_fc<=gmm_thresh[3])==1] = 1
dr_average_fc_binary[(dr_average_fc>gmm_thresh[3])] = 2

dr_average_fc_binary_index = apply(dr_average_fc_binary, 1, function(x) paste(x, collapse='_'))

pdf('trinary_index_hist.pdf')
hist((table(dr_average_fc_binary_index)), breaks=50)
dev.off()

table(dr_average_fc_binary_index)[table(dr_average_fc_binary_index)>1000]






