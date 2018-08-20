data = read.table('201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.sort.meansigsort.txt', header=F)

index = rownames(table(data[,4]))

sig_mat_all = log2(data[,-c(1:4)]+0.1)
sig_mat_all_min = min(sig_mat_all)
sig_mat_all_max = max(sig_mat_all)

euclidean_dist = c()
for (i in c(1:length(index))){
pdf(paste('kmeans.box', toString(i), '.pdf', sep=''))
sig_mat = log2(data[data[,4]==as.numeric(index[i]),-c(1:4)]+0.1)
boxplot(sig_mat, ylim=c(sig_mat_all_min, sig_mat_all_max))
lines(colMeans(sig_mat))
dev.off()
euclidean_dist[i] = sum(abs(sig_mat-colMeans(sig_mat)))
}
print('kmeans:')
print(sum(euclidean_dist))
print(sum(abs(sig_mat_all-colMeans(sig_mat_all))))
print(1-sum(euclidean_dist)/sum(abs(sig_mat_all-colMeans(sig_mat_all))))
print(head(sig_mat_all))
print(colMeans(sig_mat_all))
print(head(sig_mat_all-colMeans(sig_mat_all)))

data = read.table('201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.meansigsort.txt', header=F)

index = rownames(table(data[,4]))

sig_mat_all = log2(data[,-c(1:4)]+0.1)
sig_mat_all_min = min(sig_mat_all)
sig_mat_all_max = max(sig_mat_all)

euclidean_dist = c()
for (i in c(1:length(index))){
pdf(paste('index.box', index[i], '.pdf', sep=''))
sig_mat = log2(data[data[,4]==(index[i]),-c(1:4)]+0.1)
boxplot(sig_mat, ylim=c(sig_mat_all_min, sig_mat_all_max))
lines(colMeans(sig_mat))
dev.off()
euclidean_dist[i] = sum(abs(sig_mat-colMeans(sig_mat)))
}

print('index:')
print(sum(euclidean_dist))
print(sum(abs(sig_mat_all-colMeans(sig_mat_all))))
print(1-sum(euclidean_dist)/sum(abs(sig_mat_all-colMeans(sig_mat_all))))
