
signal_mat = read.table('201803106_AllCtCFPeaks_DiffBind.S3norm.nochrUn_random.txt', header=F)
signal_mat = as.matrix(signal_mat)

signal_mat_log2 = log2(signal_mat+0.1)
signal_mat_log2_od = signal_mat_log2
set.seed(2018)
fit <- kmeans(signal_mat_log2, 6)

kmeans_id = fit$cluster

write.table(kmeans_id, 'kmeans_id_index_vec.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

index_count = as.matrix(table(kmeans_id))

write.table(index_count, 'kmeans_id_index_count.txt', quote=FALSE, col.names=FALSE, row.names=TRUE, sep='\t')
