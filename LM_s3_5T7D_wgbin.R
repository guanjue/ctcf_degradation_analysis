#######################################
### FOR 5T7D
#######################################
print('FOR 5T7D')
all_sample_input_table = read.table('bam_bw_bedgraph/file_list.all.txt', header=F)[c(1,2,16:19, 3,4, 20:24),]
IP_list = all_sample_input_table[,1]
CTRL_list = all_sample_input_table[,2]
TP_list = apply(all_sample_input_table, 1, function(x) unlist(strsplit(x[1], split='_'))[2])
ID_list = apply(all_sample_input_table, 1, function(x) unlist(strsplit(x[1], split='_'))[1])

### get LM scale factors
LM_sf0 = read.table('s3norm_SF.5T7D.txt', header=TRUE)
LM_sf_6TP = read.table('s3norm_SF.6TP.txt', header=TRUE)

### scale 6TP & 5T7D
SF_6TP_5T7D = mean(as.matrix(LM_sf_6TP[,c(1,2,3,4)]))/mean(as.matrix(LM_sf0[,c(1,2,7,8)]))
LM_sf = LM_sf0 * SF_6TP_5T7D

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
	write.table(signal_bedgraph, paste('bin_tab/', IP_list[i], '.mm9_50.sorted.signal.fc.5T7D.bedgraph', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
}


