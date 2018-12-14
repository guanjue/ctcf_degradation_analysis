#######################################
### FOR 6 time points
#######################################
print('FOR 6 time points')
all_sample_input_table = read.table('bam_bw_bedgraph/file_list.all.txt', header=F)[c(1:15),]
IP_list = all_sample_input_table[,1]
CTRL_list = all_sample_input_table[,2]
TP_list = apply(all_sample_input_table, 1, function(x) unlist(strsplit(x[1], split='_'))[2])
ID_list = apply(all_sample_input_table, 1, function(x) unlist(strsplit(x[1], split='_'))[1])

### get LM scale factors
LM_sf = read.table('s3norm_SF.6TP.txt', header=TRUE)
### 
method_num = dim(LM_sf)[2] / length(IP_list)
sample_num = length(IP_list)

### get peak name sorted bed
peak_name_sort_pks = read.table('mm9_50.PKsorted.bed', header=FALSE)

### start
for (i in c(1:sample_num)){
	### get unadjusted signal vector
	signal_file = paste('bin_tab/', IP_list[i], '.mm9_50.sorted.signal.fc.tab', sep='')
	signal_vec = scan(signal_file)
	sf_tmp = LM_sf[5, i]
	print(sf_tmp)
	### LM adjustment
	signal_vec = signal_vec * sf_tmp
	#signal_vec = (signal_vec * sf_tmp[2] + sf_tmp[1])
	### cbind matrix
	signal_bedgraph = cbind(peak_name_sort_pks, signal_vec)
	write.table(signal_bedgraph, paste('bin_tab/', IP_list[i], '.mm9_50.sorted.signal.fc.6TP.bedgraph', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
}


