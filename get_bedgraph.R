### get parameters
args = commandArgs(trailingOnly=TRUE)

sample_list_file = args[1]
sample_list_file = 'id_list_all.txt'
sample_list = read.table(sample_list_file, header=FALSE)
sample_list = sample_list[,]

mm9_500bp_bed = read.table('mm9_500_all.sorted.bed', header=FALSE)[,c(1:3)]

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
file1 = paste(sample_list[1, 1], '.mm9_500_all.sorted.rc_per_bp.tab', sep='')
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

rep_sig_s3norm_mat_ave = c()
rep_sig_TSnorm_mat_ave = c()

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

	print(sample_list[,3])
	for (j in c(1:length(rep_list))){
		AB_used_tmp_j = AB_used_tmp[j,]
		AB_TS_used_tmp_j = AB_TS_used_tmp[j,]

		rep1 = paste(rep_list[j], '.mm9_500_all.sorted.rc_per_bp.tab', sep='')
		input = paste(input_list, '.mm9_500_all.input_sig.macs.tab', sep='')

		d_chipseq_1 = scan(rep1)*1000
		d_chipseq_1_input = scan(input)*1000
		d_chipseq_1_input[is.na(d_chipseq_1_input)] = 0

		### get FC
		d_chipseq_1_fc = (d_chipseq_1+smallnum) / (d_chipseq_1_input+smallnum)
		### s3norm
		d_chipseq_1_fc_s3norm = as.numeric(AB_used_tmp_j[1]) * d_chipseq_1_fc ^ as.numeric(AB_used_tmp_j[2])
		### TSnorm
		d_chipseq_1_fc_TSnorm = as.numeric(AB_TS_used_tmp_j[1]) * d_chipseq_1_fc

		### save to ave matrix
		d_chipseq_1_mat_tmp = cbind(d_chipseq_1_mat_tmp, d_chipseq_1)
		### save bedgraph
		d_chipseq_1_fc_s3norm_bedgraph = cbind(mm9_500bp_bed, d_chipseq_1_fc_s3norm)
		write.table(d_chipseq_1_fc_s3norm_bedgraph, paste(rep_list[j], '.mm9_500_all.sorted.rc_per_bp.tab', '.s3norm.bedgraph', sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
		d_chipseq_1_fc_TSnorm_bedgraph = cbind(mm9_500bp_bed, d_chipseq_1_fc_TSnorm)
		write.table(d_chipseq_1_fc_TSnorm_bedgraph, paste(rep_list[j], '.mm9_500_all.sorted.rc_per_bp.tab', '.TSnorm.bedgraph', sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
	}
	d_chipseq_1_ave_fc_tmp = (rowMeans(d_chipseq_1_mat_tmp)+smallnum)/(d_chipseq_1_input+smallnum)
	### ave s3
	rep_sig_s3norm_mat_ave_tmp = as.numeric(AB_used_tmp_j[1]) * d_chipseq_1_ave_fc_tmp ^ as.numeric(AB_used_tmp_j[2])
	### ave TS
	rep_sig_TSnorm_mat_ave_tmp = as.numeric(AB_TS_used_tmp_j[1]) * d_chipseq_1_ave_fc_tmp 
	### write average bedgraph
	d_chipseq_1_fc_s3norm_ave_bedgraph = cbind(mm9_500bp_bed, rep_sig_s3norm_mat_ave_tmp)
	d_chipseq_1_fc_TSnorm_ave_bedgraph = cbind(mm9_500bp_bed, rep_sig_TSnorm_mat_ave_tmp)
	write.table(d_chipseq_1_fc_s3norm_ave_bedgraph, paste(tp_type[i], '.mm9_500_all.sorted.rc_per_bp.tab', '.s3norm.ave.bedgraph', sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
	write.table(d_chipseq_1_fc_TSnorm_ave_bedgraph, paste(tp_type[i], '.mm9_500_all.sorted.rc_per_bp.tab', '.s3norm.ave.bedgraph', sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
}





time sort -k1,1 -k2,2n 1579_0A.mm9_500_all.sorted.rc_per_bp.tab.s3norm.bedgraph > 1579_0A.mm9_500_all.sorted.rc_per_bp.tab.s3norm.sort.bedgraph
time /storage/home/gzx103/group/software/ucsc/bedGraphToBigWig 1579_0A.mm9_500_all.sorted.rc_per_bp.tab.s3norm.sort.bedgraph /storage/home/gzx103/group/genome/mm9/mm9.chrom.sizes 1579_0A.mm9_500_all.sorted.rc_per_bp.tab.s3norm.bw
rm 1579_0A.mm9_500_all.sorted.rc_per_bp.tab.s3norm.sort.bedgraph



