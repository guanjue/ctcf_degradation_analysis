ave = as.matrix(read.table('rep_sig_s3norm_mat_ave.txt', header=TRUE))
rep = as.matrix(read.table('rep_sig_s3norm_mat.txt', header=TRUE))
bed = as.matrix(read.table('allmerged.sorted.bed', header=FALSE))

set.seed(2018)

used_id = sample(dim(ave)[1], 10000)
ave_r = ave[used_id,]
rep_r = rep[used_id,]
bed_r = bed[used_id,4]

mat_DPGP = cbind(bed_r, ave_r)

write.table(ave_r, 'rep_sig_s3norm_mat_ave.r.txt', quote=F, sep='\t', col.names=TRUE, row.names=F)
write.table(mat_DPGP[,-c(8:11)], 'rep_sig_s3norm_mat_ave.r.DPGP.txt', quote=F, sep='\t', col.names=TRUE, row.names=F)

write.table(ave_r[,-c(7:10)], 'rep_sig_s3norm_mat_ave.nodrb_tri.r.txt', quote=F, sep='\t', col.names=F, row.names=F)

write.table(rep_r, 'rep_sig_s3norm_mat.r.txt', quote=F, sep='\t', col.names=TRUE, row.names=F)



ave = as.matrix(read.table('rep_sig_TSnorm_mat_ave.txt', header=TRUE))
rep = as.matrix(read.table('rep_sig_TSnorm_mat.txt', header=TRUE))
bed = as.matrix(read.table('allmerged.sorted.bed', header=FALSE))

set.seed(2018)

used_id = sample(dim(ave)[1], 10000)
ave_r = ave[used_id,]
rep_r = rep[used_id,]
bed_r = bed[used_id,4]

mat_DPGP = cbind(bed_r, ave_r)

write.table(ave_r, 'rep_sig_TSnorm_mat_ave.r.txt', quote=F, sep='\t', col.names=TRUE, row.names=F)
write.table(mat_DPGP[,-c(8:11)], 'rep_sig_TSnorm_mat_ave.r.DPGP.txt', quote=F, sep='\t', col.names=TRUE, row.names=F)

write.table(ave_r[,-c(7:10)], 'rep_sig_TSnorm_mat_ave.nodrb_tri.r.txt', quote=F, sep='\t', col.names=F, row.names=F)

write.table(rep_r, 'rep_sig_TSnorm_mat.r.txt', quote=F, sep='\t', col.names=TRUE, row.names=F)


head -1 rep_sig_s3norm_mat_ave.r.DPGP.txt > rep_sig_s3norm_mat_ave.percent.r.DPGP.txt
tail -n+2 rep_sig_s3norm_mat_ave.r.DPGP.txt | awk -F '\t' -v OFS='\t' '{print $1, $2*100/$2, $3*100/$2, $4*100/$2, $5*100/$2, $6*100/$2, $7*100/$2}' >> rep_sig_s3norm_mat_ave.percent.r.DPGP.txt


head -1 rep_sig_TSnorm_mat_ave.r.DPGP.txt > rep_sig_TSnorm_mat_ave.percent.r.DPGP.txt
tail -n+2 rep_sig_TSnorm_mat_ave.r.DPGP.txt | awk -F '\t' -v OFS='\t' '{print $1, ($2+1)*100/($2+1), ($3+1)*100/($2+1), ($4+1)*100/($2+1), ($5+1)*100/($2+1), ($6+1)*100/($2+1), ($7+1)*100/($2+1)}' >> rep_sig_TSnorm_mat_ave.percent.r.DPGP.txt

head -1 rep_sig_TSnorm_mat_ave.r.DPGP.txt > rep_sig_TSnorm_mat_ave.percent.all.DPGP.txt
tail -n+2 rep_sig_TSnorm_mat_ave.txt | awk -F '\t' -v OFS='\t' '{print ($1+1)*100/($1+1), ($2+1)*100/($1+1), ($3+1)*100/($1+1), ($4+1)*100/($1+1), ($5+1)*100/($1+1), ($6+1)*100/($1+1)}' > rep_sig_TSnorm_mat_ave.percent.r.DPGP.txt.tmp
paste allmerged.sorted.bed rep_sig_TSnorm_mat_ave.percent.r.DPGP.txt.tmp | awk -F '\t' -v OFS='\t' '{print $4,$5,$6,$7,$8,$9,$10}' >> rep_sig_TSnorm_mat_ave.percent.all.DPGP.txt

