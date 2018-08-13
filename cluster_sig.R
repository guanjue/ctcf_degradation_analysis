cp IDEAS_run_linear/run_IDEAS_result/run_IDEAS.wg.sort.state.bed ./

tail -n+2 run_IDEAS.wg.sort.state.bed | awk -F '\t' -v OFS='\t' '{if ($4!="0_0_0_0") print $0}' | sort -k1,1 -k2,2n > run_IDEAS.wg.sort.state.pk.bed


bedtools intersect -a 201803106_AllCtCFPeaks_DiffBind.bed -b run_IDEAS.wg.sort.state.pk.bed -wa -u > diffbind.ideas.pk.bed

bedtools merge -i  run_IDEAS.wg.sort.state.pk.bed >  run_IDEAS.wg.sort.state.pk.merge.bed



.2r_nbp.fisher_p.200bp.txt


paste mm9_200.sort.bed $ct'.2r_nbp.fisher_p.200bp.txt' > $ct'.2r_nbp.fisher_p.200bp.bedgraph'

time /storage/home/gzx103/group/software/ucsc/bedGraphToBigWig $ct'.2r_nbp.fisher_p.200bp.normed.txt' /storage/home/gzx103/group/projects/vision/input_norm/mm9.chrom.sizes $ct'.S3norm.bw'



paste WT.2r_nbp.fisher_p.200bp.normed.txt 0hr.2r_nbp.fisher_p.200bp.normed.txt 4hr.2r_nbp.fisher_p.200bp.normed.txt 6hr.2r_nbp.fisher_p.200bp.normed.txt > mm9.200bp.signal.mat.txt

signal_mat = read.table('mm9.200bp.signal.mat.txt', header=F)

signal_mat = as.matrix(signal_mat)

signal_mat_log2 = log2(signal_mat+0.1)

png('hist.S3norm.png')
plot(density(signal_mat[signal_mat!=0]))
dev.off()

png('hist.S3norm.log2.png')
hist(signal_mat_log2, breaks = 50)
dev.off()

set.seed(2018)
sample_id = sample(dim(signal_mat_log2)[1], 100000)
signal_mat_log2_s = signal_mat_log2[sample_id]

library(mixtools)
gmm_k = 5
mixmdl = normalmixEM(signal_mat_log2_s, k = gmm_k)
png('tp.sig.log2.gmm.png')
plot(mixmdl, which = 2)
dev.off()

mixmdl$mu


multmixmodel.sel(signal_mat_log2_s)

signal_mat_log2_s = as.matrix(signal_mat_log2_s)

set.seed(2018)
library(mixtools)
gmm_k = 4
mixmdl = normalmixEM(signal_mat_log2_s[signal_mat_log2_s>2], k = gmm_k)
png('tp.sig.log2.gmm.png')
plot(mixmdl, which = 2)
lines(density(signal_mat_log2_s[signal_mat_log2_s>2]), lty = 2, lwd = 2)
dev.off()


c = rowSums(signal_mat)

signal_mat_label = signal_mat

signal_mat_label = signal_mat_label[c]
signal_mat_label[signal_mat<=2]=0
signal_mat_label[(signal_mat>2)*(signal_mat<=2.6)==1]=0
signal_mat_label[(signal_mat>2.6)*(signal_mat<=3.8)==1]=1
signal_mat_label[(signal_mat>3.8)*(signal_mat<=6.3)==1]=2
signal_mat_label[(signal_mat>6.3)]=3

signal_mat_label_fast = signal_mat_label[sample_id,]
signal_mat_label_vec = apply(signal_mat_label_fast,1,function(x) paste(x[1], x[2], x[3], x[4], sep='_'))

tl = as.matrix(table(signal_mat_label_vec))

tl[order(rownames(tl))]



