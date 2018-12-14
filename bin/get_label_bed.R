#data = as.data.frame(read.table('mm9_200.sort.bed', header=FALSE))
#data = as.data.frame(read.table('mm9_200.sort.mini.bed', header=FALSE))
data = as.data.frame(read.table('201803106_AllCtCFPeaks_DiffBind.bed', header=FALSE))

data_label = c(1:dim(data)[1])

data_new = cbind(data, data_label)

#write.table(data_new, 'mm9_200.sort.label.bed', sep=' ', quote=FALSE, col.names=FALSE, row.names=FALSE)
#write.table(data_new, 'mm9_200.sort.mini.label.bed', sep=' ', quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(data_new, '201803106_AllCtCFPeaks_DiffBind.label.bed', sep=' ', quote=FALSE, col.names=FALSE, row.names=FALSE)
