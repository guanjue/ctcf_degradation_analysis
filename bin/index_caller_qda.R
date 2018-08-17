library(MASS)

### read signal matrix
signal_mat = read.table('201803106_AllCtCFPeaks_DiffBind.S3norm.nochrUn_random.txt', header=F)
signal_mat = as.data.frame(signal_mat)
signal_mat_log2 = log2(signal_mat+0.1)
### convert to fold change
signal_mat_log2_fc = signal_mat_log2-(signal_mat_log2[,2])

### read all pk index vec
signal_mat_index_od = as.matrix(read.table('signal_mat_index_vec.txt', header=F))
signal_mat_index = signal_mat_index_od

### read all index
index_count_table = as.matrix(read.table('index_count_all.txt', header=F))
index_count_thresh = 200
index_count_table_small = index_count_table[as.numeric(index_count_table[,2])<=index_count_thresh,1]

for (small_index in index_count_table_small){
	print(small_index)
	signal_mat_index[signal_mat_index_od==toString(small_index)] = 'X_X_X_X'
}


### as matrix
signal_mat_index_sig = cbind(signal_mat_index, signal_mat_log2[,2], signal_mat_log2_fc[,c(1,3,4)])
colnames(signal_mat_index_sig) = c('index', '0hr_sig', 'WT_fc', '4hr_fc', '6hr_fc')

index_vec = signal_mat_index_sig[,1]


### train qda
print(as.matrix(table(index_vec)))

for (i in c(1:1000)){
	index_count_old = as.matrix(table(index_vec))
	qda_model = qda(signal_mat_index_sig[,c(2:5)], index_vec)
	### predict based on trained qda model
	predicted_index = predict(qda_model, signal_mat_index_sig[,c(2:5)])$class
	index_vec = predicted_index
	index_count_new = as.matrix(table(index_vec))
	dif = sum(abs(index_count_new-index_count_old))
	if (i%%10==0){
		print(paste('iter: ', toString(i), sep=''))
		print(dif)
	}
	if (dif==0){
		print(paste('iter: ', toString(i), sep=''))
		print(dif)
		break
	}
}

### merge all small index set
predicted_index = as.matrix(predicted_index)

new_index_count_table_small = rownames(index_count_new)[as.numeric(index_count_new)<index_count_thresh]
predicted_index_output = predicted_index
for (small_index in new_index_count_table_small){
	print(small_index)
	print(sum(predicted_index==toString(small_index)))
	predicted_index_output[predicted_index==toString(small_index)] = 'X_X_X_X'
}

index_count_new_nosmall = as.matrix(table(predicted_index_output))

print(index_count_new_nosmall)

write.table(predicted_index_output, 'signal_mat_index_vec_ic.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
write.table(index_count_new_nosmall, 'index_count_enrich_ic.txt', quote=FALSE, col.names=FALSE, row.names=TRUE, sep='\t')




