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
AB_1 = read.table('rep1_s3norm_AB.txt', sep='\t', header=TRUE)
AB_2 = read.table('rep2_s3norm_AB.txt', sep='\t', header=TRUE)
AB_average = read.table('average_s3norm_AB.txt', sep='\t', header=TRUE)

i = 1

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
	AB_1_i = AB_1[i,]
	AB_2_i = AB_1[i,]
	AB_average_i = AB_average[i,]
	### normalize
	rep1_sig_fc = rep1_sig/input_sig
	rep2_sig_fc = rep2_sig/input_sig
	rep1_sig_s3norm = (rep1_sig_fc)^AB_1_i[1,2] * AB_1_i[1,1]
	rep2_sig_s3norm = (rep2_sig_fc)^AB_2_i[1,2] * AB_2_i[1,1]
	rep1_sig_s3norm_mat = cbind(rep1_sig_s3norm_mat, rep1_sig_s3norm)
	rep2_sig_s3norm_mat = cbind(rep2_sig_s3norm_mat, rep2_sig_s3norm)
	### qtpcr
	rep1_sig_fc_qtpcr = rep1_sig_qtpcr/input_sig_qtpcr
	rep2_sig_fc_qtpcr = rep2_sig_qtpcr/input_sig_qtpcr
	rep1_sig_s3norm_qtpcr = (rep1_sig_fc_qtpcr)^AB_average_i[1,2] * AB_average_i[1,1]
	rep2_sig_s3norm_qtpcr = (rep2_sig_fc_qtpcr)^AB_average_i[1,2] * AB_average_i[1,1]
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



