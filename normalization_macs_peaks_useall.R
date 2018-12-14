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

smallnum = 1
### normalization
for (i in c(1:length(tp_type))){
	rep_list = sample_list[sample_list[,3] == tp_type[i],1]
	input_list = sample_list[sample_list[,3] == tp_type[i],2][1]
	AB_used_tmp = AB_used[AB_used[,3]==tp_type[i],]
	AB_average_used_tmp = AB_average_used[AB_average_used[,3]==tp_type[i],]
	AB_TS_used_tmp = AB_TS_used[AB_TS_used[,3]==tp_type[i],]
	AB_TS_average_used_tmp = AB_TS_average_used[AB_TS_average_used[,3]==tp_type[i],]

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
	}
}

colnames(rep_sig_s3norm_mat) = sample_list[,3]
colnames(rep_sig_TSnorm_mat) = sample_list[,3]
colnames(rep_sig_s3norm_mat_nonpk) = sample_list[,3]
colnames(rep_sig_TSnorm_mat_nonpk) = sample_list[,3]

rep_sig_s3norm_mat_ave = c()
rep_sig_TSnorm_mat_ave = c()
rep_sig_s3norm_mat_nonpk_ave = c()
rep_sig_TSnorm_mat_nonpk_ave = c()
for (i in c(1:length(tp_type))){
	rep_sig_s3norm_mat_ave_tmp = rowMeans(rep_sig_s3norm_mat[, colnames(rep_sig_s3norm_mat)==tp_type[i]])
	rep_sig_s3norm_mat_ave = cbind(rep_sig_s3norm_mat_ave, rep_sig_s3norm_mat_ave_tmp)
	rep_sig_TSnorm_mat_ave_tmp = rowMeans(rep_sig_TSnorm_mat[, colnames(rep_sig_TSnorm_mat)==tp_type[i]])
	rep_sig_TSnorm_mat_ave = cbind(rep_sig_TSnorm_mat_ave, rep_sig_TSnorm_mat_ave_tmp)
	rep_sig_s3norm_mat_nonpk_ave_tmp = rowMeans(rep_sig_s3norm_mat_nonpk[, colnames(rep_sig_s3norm_mat_nonpk)==tp_type[i]])
	rep_sig_s3norm_mat_nonpk_ave = cbind(rep_sig_s3norm_mat_nonpk_ave, rep_sig_s3norm_mat_nonpk_ave_tmp)
	rep_sig_TSnorm_mat_nonpk_ave_tmp = rowMeans(rep_sig_TSnorm_mat_nonpk[, colnames(rep_sig_TSnorm_mat_nonpk)==tp_type[i]])
	rep_sig_TSnorm_mat_nonpk_ave = cbind(rep_sig_TSnorm_mat_nonpk_ave, rep_sig_TSnorm_mat_nonpk_ave_tmp)
}

colnames(rep_sig_s3norm_mat_ave) = tp_type
colnames(rep_sig_TSnorm_mat_ave) = tp_type
colnames(rep_sig_s3norm_mat_nonpk_ave) = tp_type
colnames(rep_sig_TSnorm_mat_nonpk_ave) = tp_type


write.table(rep_sig_s3norm_mat, 'rep_sig_s3norm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_s3norm_mat_nonpk, 'rep_sig_s3norm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_s3norm_mat_ave, 'rep_sig_s3norm_mat_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_s3norm_mat_nonpk_ave, 'rep_sig_s3norm_mat_nopk_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(rep_sig_TSnorm_mat, 'rep_sig_TSnorm_mat.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_TSnorm_mat_nonpk, 'rep_sig_TSnorm_mat_nopk.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_TSnorm_mat_ave, 'rep_sig_TSnorm_mat_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(rep_sig_TSnorm_mat_nonpk_ave, 'rep_sig_TSnorm_mat_nopk_ave.txt', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)



