#paste WT.qt_pcr_sig.txt 0hr.qt_pcr_sig.txt 4hr.qt_pcr_sig.txt 6hr.qt_pcr_sig.txt >  ALL.qt_pcr_sig.txt
#paste WT.2r_nbp.fisher_p.qtpcr.txt 0hr.2r_nbp.fisher_p.qtpcr.txt 4hr.2r_nbp.fisher_p.qtpcr.txt 6hr.2r_nbp.fisher_p.qtpcr.txt >  ALL.chipseq_sig.txt

#paste WT.nbp_2r_bgadj.rep1.qtpcr.txt 0hr.nbp_2r_bgadj.rep1.qtpcr.txt 4hr.nbp_2r_bgadj.rep1.qtpcr.txt 6hr.nbp_2r_bgadj.rep1.qtpcr.txt >  ALL.chipseq_sig.repsep.txt
#paste WT.nbp_2r_bgadj.rep2.qtpcr.txt 0hr.nbp_2r_bgadj.rep2.qtpcr.txt 4hr.nbp_2r_bgadj.rep2.qtpcr.txt 6hr.nbp_2r_bgadj.rep2.qtpcr.txt >>  ALL.chipseq_sig.repsep.txt
#paste WT.qt_pcr_sig.1.txt 0hr.qt_pcr_sig.1.txt 4hr.qt_pcr_sig.1.txt 6hr.qt_pcr_sig.1.txt >  ALL.qt_pcr_sig.repsep.txt
#paste WT.qt_pcr_sig.2.txt 0hr.qt_pcr_sig.2.txt 4hr.qt_pcr_sig.2.txt 6hr.qt_pcr_sig.2.txt >>  ALL.qt_pcr_sig.repsep.txt


### read signal
qtpcr=as.matrix(read.table('ALL.qt_pcr_sig.txt', header=FALSE))
chip_fihserp=as.matrix(read.table('ALL.chipseq_sig.txt', header=FALSE))

small_num = 1e-3
upperlim = 323

### check distribution
png('hist.qtpcr.png')
plot(density(qtpcr))
dev.off()
png('hist.chipseq.png')
plot(density(chip_fihserp))
dev.off()

### check 2d relationship
png('test1.nolog2chipseq.png')
plot(as.vector(qtpcr), as.vector((chip_fihserp)))
dev.off()
png('test1.log2chipseq.png')
plot(as.vector(qtpcr), as.vector(log2(chip_fihserp)))
dev.off()
png('test2.log2qtpcr.png')
plot(as.vector(log2(qtpcr+small_num)), as.vector((chip_fihserp)))
dev.off()
png('test3.log2both.png')
pr2 = cor(as.vector(log2(chip_fihserp)), as.vector(log2(qtpcr+small_num)))
print(pr2)
plot(as.vector(log2(qtpcr+small_num)), as.vector(log2(chip_fihserp)), main=toString(pr2))
lim_fit = lm(as.vector(log2(chip_fihserp))~as.vector(log2(qtpcr+small_num)))
lines(as.vector(log2(qtpcr+small_num)) ,lim_fit$fitted.values)
dev.off()

### normalize qtPCR
log2qtpcr = log2(qtpcr+small_num)
log2chipseq = log2(chip_fihserp+small_num)
log2qtpcr_norm = ( log2qtpcr - mean(log2qtpcr) ) / sd(log2qtpcr) * sd(log2chipseq) + mean(log2chipseq)
#log2qtpcr_norm = log2qtpcr * lim_fit$coefficients[2] + lim_fit$coefficients[1]

### plot scatterplot after normalization
png('normqtpcr_vs_chipseq.log2both.png')
qtPCR_region_cor = cor(log2qtpcr_norm, log2chipseq)
plot(log2qtpcr_norm, log2chipseq, xlim=c(-3, 9), ylim=c(-3, 9), main=toString(c(round(qtPCR_region_cor[1,1],3),round(qtPCR_region_cor[2,2],3),round(qtPCR_region_cor[3,3],3),round(qtPCR_region_cor[4,4],3))))
abline(0,1,col='red')
dev.off()


for (i in c(1:4)){
	png(paste('normqtpcr_vs_chipseq.log2both', toString(i), '.png'))
	plot(log2qtpcr_norm[,i], log2chipseq[,i], xlim=c(-3, 9), ylim=c(-3, 9))
	abline(0,1,col='red')
	dev.off()
}

### read DNA region type
regions = read.table('qt_pcr_sig.sort.txt', header=TRUE)
regions_name = apply(regions, 1, function(x) toString(x[1]))
regions_name2 = c(regions_name,regions_name)

### read chip-seq signal (whole genome)
chipseq_od_WT = scan('WT.2r_nbp.fisher_p.200bp.txt')
chipseq_od_0hr = scan('0hr.2r_nbp.fisher_p.200bp.txt')
chipseq_od_4hr = scan('4hr.2r_nbp.fisher_p.200bp.txt')
chipseq_od_6hr = scan('6hr.2r_nbp.fisher_p.200bp.txt')
chipseq_od_WT_pk = p.adjust(10^(-chipseq_od_WT), 'fdr')
chipseq_od_0hr_pk = p.adjust(10^(-chipseq_od_0hr), 'fdr')
chipseq_od_4hr_pk = p.adjust(10^(-chipseq_od_4hr), 'fdr')
chipseq_od_6hr_pk = p.adjust(10^(-chipseq_od_6hr), 'fdr')
chipseq_od_mat_all = cbind(chipseq_od_WT, chipseq_od_0hr, chipseq_od_4hr, chipseq_od_6hr)
chipseq_od_mat_all_pk = cbind(chipseq_od_WT_pk, chipseq_od_0hr_pk, chipseq_od_4hr_pk, chipseq_od_6hr_pk)

### set random id
set.seed(2018)
sample_id = sample(length(chipseq_od_WT), 1000000)
sample_id_loess = sample(length(sample_id), 50000)

### plot qtPCR vs chipseq (whole genome)
for (i in c(1:4)){
print(i)
png(paste('test.without_norm.chipseq_qtpcr.', toString(i), '.png', sep=''))
chip_tmp1 = log2(chipseq_od_mat_all[sample_id,2]+small_num)
chip_tmp2 = log2(chipseq_od_mat_all[sample_id,i]+small_num)
chip_tmp1_qtpcr_region = log2chipseq[,2]
chip_tmp2_qtpcr_region = log2chipseq[,i]
chip_tmp1_pk = chipseq_od_mat_all_pk[sample_id,2] < 0.05
chip_tmp2_pk = chipseq_od_mat_all_pk[sample_id,i] < 0.05
cpk = chip_tmp1_pk*chip_tmp2_pk > 0
cbg = chip_tmp1_pk+chip_tmp2_pk == 0
log2qtpcr_norm_tmp1 = log2qtpcr_norm[,2]
log2qtpcr_norm_tmp2 = log2qtpcr_norm[,i]
#stable_pr2 = cor(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic'])
stable_pr2 = 1 - sum((log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']-chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'])^2) / sum((log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']-mean(log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']))^2) 
plot(chip_tmp1, chip_tmp2, xlim=c(-3, 9), ylim=c(-3, 9), col='gray', main=toString(round(stable_pr2, 3)))
#points(chip_tmp1[cpk], chip_tmp2[cpk], col='orangered')
#plot(log2qtpcr_norm_tmp1, log2qtpcr_norm_tmp2, xlim=c(-3, 9), ylim=c(-3, 9))
### plot chip-seq signal at qtPCR regions
points(chip_tmp1_qtpcr_region[regions_name=='Strong-Dynamic'], chip_tmp2_qtpcr_region[regions_name=='Strong-Dynamic'], col='blue')
points(chip_tmp1_qtpcr_region[regions_name=='Weak'], chip_tmp2_qtpcr_region[regions_name=='Weak'], col='green')
points(chip_tmp1_qtpcr_region[regions_name=='Strong-Stable'], chip_tmp2_qtpcr_region[regions_name=='Strong-Stable'], col='red')
lm_stable_chip = lm(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic']~chip_tmp1_qtpcr_region[regions_name!='Strong-Dynamic'])
abline(lm_stable_chip, col="blue")
### plot qtPCR signal at qtPCR regions
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name=='Strong-Dynamic'], col='blue', pch=16)
points(log2qtpcr_norm_tmp1[regions_name=='Weak'], log2qtpcr_norm_tmp2[regions_name=='Weak'], col='green', pch=16)
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable'], log2qtpcr_norm_tmp2[regions_name=='Strong-Stable'], col='red', pch=16)
lm_stable_qtpcr = lm(log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']~log2qtpcr_norm_tmp1[regions_name!='Strong-Dynamic'])
abline(lm_stable_qtpcr, col="red")
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable']), mean(log2qtpcr_norm_tmp2[regions_name=='Strong-Stable']), col='red', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Weak']), mean(log2qtpcr_norm_tmp2[regions_name=='Weak']), col='green', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic']), mean(log2qtpcr_norm_tmp2[regions_name=='Strong-Dynamic']), col='blue', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name!='Strong-Dynamic']), mean(log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']), col='deepskyblue', pch=16)
points(mean(chip_tmp1), mean(chip_tmp2), col='deepskyblue')
points(mean(chip_tmp1[cpk]), mean(chip_tmp2[cpk]), col='orangered')
points(mean(chip_tmp1[cbg]), mean(chip_tmp2[cbg]), col='black')
abline(0,1,col='black')
dev.off()
}



################################################################################################
### NewtonRaphsonMethod
NR_method = function(sig1_pk,sig1_bg, sig2_pk,sig2_bg, upperlim, A,B, moment, converge_thresh, numIterations){
	sig1_pk_mean = mean(sig1_pk^moment)
	sig1_bg_mean = mean(sig1_bg^moment)
	for ( i in c(1:numIterations)){
		sig2_pk_transformed = sig2_pk^(moment*B)
		sig2_bg_transformed = sig2_bg^(moment*B)
		fb = sig1_bg_mean * mean(sig2_pk_transformed) - sig1_pk_mean * mean(sig2_bg_transformed)
		dfb = moment * sig1_bg_mean * mean(log(sig2_pk) * sig2_pk_transformed) - moment * sig1_pk_mean * mean(log(sig2_bg) * sig2_bg_transformed)
		### next step
		B = B - fb / dfb
		sig2_bg_transformed = sig2_bg^(moment*B)
		A = sig1_bg_mean / mean(sig2_bg_transformed)
		print(c(i, dfb))
		print(c(A,B))
		print(abs(fb / dfb))
		print(abs(fb / dfb) < converge_thresh)
		last_AB = c(A,B)
		if (abs(fb / dfb) < converge_thresh){
			print('converged!')
			used_AB = c(A, B)
			break
		}
	}
	if (abs(fb / dfb) >= converge_thresh){
		print('NOT converged...')
		used_AB = last_AB
	}
	print('used: ')
	print(used_AB)
	return(used_AB)
}

################################################################################################
### NewtonRaphsonMethod
log2mean_method = function(sig1_pk,sig1_bg, sig2_pk,sig2_bg){
	sig1_pk_mean = mean(sig1_pk)
	sig1_bg_mean = mean(sig1_bg)
	sig2_pk_mean = mean(sig2_pk)
	sig2_bg_mean = mean(sig2_bg)
	B = (sig1_pk_mean - sig1_bg_mean) / (sig2_pk_mean - sig2_bg_mean)
	A = 2^(sig1_pk_mean - B * sig2_pk_mean)
	used_AB = c(A, B)
	print('used: ')
	print(used_AB)
	return(used_AB)
}

################################################################################################
### NewtonRaphsonMethod
Euclidean_distance = function(sig1,sig2){
	euclidean_dist = (sum((sig1 - sig2)^2))^0.5
	return(euclidean_dist)
}


test_AB_mat = c()
test_AB_mat_all = c()
Euclidean_distance_vec = c()

for (i in c(1:4)){
print(i)
chip_tmp1_pk = chipseq_od_mat_all_pk[sample_id,2] < 0.05
chip_tmp2_pk = chipseq_od_mat_all_pk[sample_id,i] < 0.05
cpk = chip_tmp1_pk*chip_tmp2_pk > 0
cbg = chip_tmp1_pk+chip_tmp2_pk == 0

chip_tmp1_bg = log2((chipseq_od_mat_all[cbg,i]+small_num))
chip_tmp2_bg = log2((chipseq_od_mat_all[cbg,2]+small_num))
qtpcr_norm_ref0_pk = log2qtpcr_norm[,2][regions_name!='Strong-Dynamic1']
qtpcr_norm_ref0_pk_all = log2qtpcr_norm[,2]
qtpcr_norm_ref_pk = log2qtpcr_norm[,i][regions_name!='Strong-Dynamic1']
qtpcr_norm_ref_pk_all = log2qtpcr_norm[,i]
#qtpcr_norm_ref_pk = (qtpcr_norm_ref_pk-mean(qtpcr_norm_ref_pk))/sd(qtpcr_norm_ref_pk) * sd(qtpcr_norm_ref0_pk) + mean(qtpcr_norm_ref0_pk)
qtpcr_norm_ref_bg = mean(chip_tmp2_bg)
chip_fihserp_tar_pk = log2chipseq[,i][regions_name!='Strong-Dynamic1']
chip_fihserp_tar_pk_all = log2chipseq[,i]
chip_fihserp_tar_bg = mean(chip_tmp1_bg)
test_AB = log2mean_method(qtpcr_norm_ref_pk,qtpcr_norm_ref_bg, chip_fihserp_tar_pk,chip_fihserp_tar_bg)
test_AB_all = log2mean_method(qtpcr_norm_ref_pk_all,qtpcr_norm_ref_bg, chip_fihserp_tar_pk_all,chip_fihserp_tar_bg)

test_AB_mat = cbind(test_AB_mat, test_AB)
test_AB_mat_all = cbind(test_AB_mat_all, test_AB_all)
Euclidean_distance_vec[i] = (sum((qtpcr_norm_ref_pk - chip_fihserp_tar_pk_all[regions_name!='Strong-Dynamic'])^2))^0.5
print('correlation: ')
print(cor(qtpcr_norm_ref_pk, chip_fihserp_tar_pk_all[regions_name!='Strong-Dynamic1']))
print('Euclidean distance: ')
print((sum((qtpcr_norm_ref_pk - chip_fihserp_tar_pk_all[regions_name!='Strong-Dynamic1'])^2))^0.5)
}



### plot qtPCR vs chipseq normed (whole genome)
timepoints = c('WT', '0hr', '4hr', '6hr')
for (i in c(1:4)){
print(i)
#i=1
A1 = test_AB_mat[1,2]
B1 = test_AB_mat[2,2]
A2 = test_AB_mat[1,i]
B2 = test_AB_mat[2,i]
A1_useall = test_AB_mat_all[1,2]
B1_useall = test_AB_mat_all[2,2]
A2_useall = test_AB_mat_all[1,i]
B2_useall = test_AB_mat_all[2,i]
###
png(paste('test.norm.chipseq_qtpcr.s3norm.compare0HR.', toString(i), '.png', sep=''))
###
chip_tmp1 = log2(chipseq_od_mat_all[sample_id,2]+small_num)
chip_tmp2 = log2(chipseq_od_mat_all[sample_id,i]+small_num)
chip_tmp1_qtpcr_region = log2chipseq[,2]
chip_tmp2_qtpcr_region = log2chipseq[,i]
chip_tmp1_pk = chipseq_od_mat_all_pk[sample_id,2] < 0.05
chip_tmp2_pk = chipseq_od_mat_all_pk[sample_id,i] < 0.05
cpk = chip_tmp1_pk*chip_tmp2_pk > 0
cbg = chip_tmp1_pk+chip_tmp2_pk == 0
###
chip_tmp1 = log2(A1 * (chipseq_od_mat_all[sample_id,2]+small_num)^B1)
chip_tmp1_useall = log2(A1_useall * (chipseq_od_mat_all[sample_id,2]+small_num)^B1_useall)
chipseq_od_mat_all_tmp = chipseq_od_mat_all[,i]
chip_tmp2 = log2(A2 * (chipseq_od_mat_all_tmp[sample_id]+small_num)^B2)
chip_tmp2_useall = log2(A2_useall * (chipseq_od_mat_all_tmp[sample_id]+small_num)^B2_useall)
chip_tmp2_all =  (A2 * (chipseq_od_mat_all_tmp+small_num)^B2)
chip_tmp2_all[chip_tmp2_all>upperlim] = upperlim
chip_tmp2_all[chipseq_od_mat_all_tmp==0] = 0.0
#chip_tmp2_all_log2 = log2(chip_tmp2_all)
write.table(chip_tmp2_all, paste(timepoints[i], '.2r_nbp.fisher_p.200bp.normed.txt', sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
#write.table(chip_tmp2_all_log2, paste(timepoints[i+1], '.2r_nbp.fisher_p.200bp.normed.log2.txt', sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
###
chip_tmp1_qtpcr_region = log2(A1 * (chip_fihserp[,2]+small_num)^B1)
chip_tmp2_qtpcr_region =  log2(A2 * (chip_fihserp[,i]+small_num)^B2)
chip_tmp1_qtpcr_region_useall = log2(A1_useall * (chip_fihserp[,2]+small_num)^B1_useall)
chip_tmp2_qtpcr_region_useall =  log2(A2_useall * (chip_fihserp[,i]+small_num)^B2_useall)
chip_tmp2_qtpcr_region_od =  (chip_fihserp[,i]+small_num)
###
chip_tmp1_pk = chipseq_od_mat_all_pk[sample_id,2] < 0.05
chip_tmp2_pk = chipseq_od_mat_all_pk[sample_id,i] < 0.05
cpk = chip_tmp1_pk*chip_tmp2_pk > 0
cbg = chip_tmp1_pk+chip_tmp2_pk == 0
plot(chip_tmp1, chip_tmp2, xlim=c(-3, 9), ylim=c(-3, 9), col='gray')
data_tmp = as.data.frame(cbind(chip_tmp1, chip_tmp2)[sample_id_loess,])
data_tmp = data_tmp[order(data_tmp[,1]),]
colnames(data_tmp) = c('x', 'y')
loess_fit = loess(y ~ x, data_tmp, span=0.2, degree=2)
lines(data_tmp$x, predict(loess_fit), col = "blue", lty=2)
log2qtpcr_norm_tmp1_stable = log2qtpcr_norm[regions_name!='Strong-Dynamic',2]
log2qtpcr_norm_tmp2_stable = log2qtpcr_norm[regions_name!='Strong-Dynamic',i]
log2qtpcr_norm_tmp1 = log2qtpcr_norm[,2]
log2qtpcr_norm_tmp2 = log2qtpcr_norm[,i]
chip_tmp2_qtpcr_region_scale = log2(chip_tmp2_qtpcr_region_od) - mean((chip_tmp2_qtpcr_region_od[regions_name!='Strong-Dynamic'])) + mean(2^log2qtpcr_norm_tmp2)
#log2qtpcr_norm_tmp2 = (log2qtpcr_norm[,i+1] - mean(log2qtpcr_norm_tmp2_stable))/sd(log2qtpcr_norm_tmp2_stable) * sd(log2qtpcr_norm_tmp1_stable)+mean(log2qtpcr_norm_tmp1_stable)
#plot(log2qtpcr_norm_tmp1, log2qtpcr_norm_tmp2, xlim=c(-3, 9), ylim=c(-3, 9))
points(chip_tmp1_qtpcr_region[regions_name=='Strong-Dynamic'], chip_tmp2_qtpcr_region[regions_name=='Strong-Dynamic'], col='blue')
points(chip_tmp1_qtpcr_region[regions_name=='Weak'], chip_tmp2_qtpcr_region[regions_name=='Weak'], col='green')
points(chip_tmp1_qtpcr_region[regions_name=='Strong-Stable'], chip_tmp2_qtpcr_region[regions_name=='Strong-Stable'], col='red')
###
#points(mean(chip_tmp1_qtpcr_region[regions_name!='Strong-Dynamic']), mean(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic']), col='purple', pch=16)
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name=='Strong-Dynamic'], col='blue', pch=16)
points(log2qtpcr_norm_tmp1[regions_name=='Weak'], log2qtpcr_norm_tmp2[regions_name=='Weak'], col='green', pch=16)
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable'], log2qtpcr_norm_tmp2[regions_name=='Strong-Stable'], col='red', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable']), mean(log2qtpcr_norm_tmp2[regions_name=='Strong-Stable']), col='red', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Weak']), mean(log2qtpcr_norm_tmp2[regions_name=='Weak']), col='green', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic']), mean(log2qtpcr_norm_tmp2[regions_name=='Strong-Dynamic']), col='blue', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name!='Strong-Dynamic']), mean(log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']), col='deepskyblue', pch=16)
points(mean(chip_tmp1), mean(chip_tmp2), col='black')
#points(mean(chip_tmp1[cpk]), mean(chip_tmp2[cpk]), col='yellow', pch=16)
#points(mean(chip_tmp1[cbg]), mean(chip_tmp2[cbg]), col='black', pch=16)
abline(0,1,col='black')
dev.off()
print('correlation: ')
print(cor(chip_tmp2_qtpcr_region, log2qtpcr_norm_tmp2))
print(cor(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']))
print(cor(chip_tmp2_qtpcr_region_useall, log2qtpcr_norm_tmp2))
print(cor(chip_tmp2_qtpcr_region_useall[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']))
print(cor(chip_tmp2_qtpcr_region_od, log2qtpcr_norm_tmp2))
print(cor(chip_tmp2_qtpcr_region_od[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']))
print(cor(chip_tmp2_qtpcr_region_scale, log2qtpcr_norm_tmp2))
print(cor(chip_tmp2_qtpcr_region_scale[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']))
print(Euclidean_distance(chip_tmp2_qtpcr_region, log2qtpcr_norm_tmp2))
print(Euclidean_distance(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']))
print(Euclidean_distance(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic'])/Euclidean_distance_vec[i])
print(Euclidean_distance(chip_tmp2_qtpcr_region_scale[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic'])/Euclidean_distance_vec[i])

###
}



print(test_AB_mat)
rownames(test_AB_mat) = c('A', 'B')
colnames(test_AB_mat) = c('WT', '0hr', '4hr', '6hr')
write.table(test_AB_mat, 'test_AB_mat.txt', quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)









png(paste('test.norm.chipseq_qtpcr.qtpcr_region.s3norm.compare0HR.', toString(i), '.png', sep=''))
chip_tmp1 = log2(test_AB[1] * (chip_fihserp[,i]+small_num)^test_AB[2])
chip_tmp1_all = log2(test_AB[1] * (chipseq_od_mat_all[sample_id,i]+small_num)^test_AB[2] )
log2qtpcr_norm_tmp1 = log2qtpcr_norm[,i]
plot(log2qtpcr_norm_tmp1, chip_tmp1, xlim=c(-3, 9), ylim=c(-3, 9), col='gray')
#plot(log2qtpcr_norm_tmp1, log2qtpcr_norm_tmp2, xlim=c(-3, 9), ylim=c(-3, 9))
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic'], chip_tmp1[regions_name=='Strong-Dynamic'], col='blue')
points(log2qtpcr_norm_tmp1[regions_name=='Weak'], chip_tmp1[regions_name=='Weak'], col='green')
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable'], chip_tmp1[regions_name=='Strong-Stable'], col='red')
points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable']), mean(chip_tmp1[regions_name=='Strong-Stable']), col='red', pch=16)
points(mean(log2qtpcr_norm_tmp1[regions_name=='Weak']), mean(chip_tmp1[regions_name=='Weak']), col='green', pch=16)
points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic']), mean(chip_tmp1[regions_name=='Strong-Dynamic']), col='blue', pch=16)
points(mean(log2qtpcr_norm_tmp1[regions_name!='Strong-Dynamic']), mean(chip_tmp1[regions_name!='Strong-Dynamic']), col='deepskyblue', pch=16)
points(mean(log2qtpcr_norm_tmp1), mean(chip_tmp1), col='black', pch=16)
points(mean(chip_tmp2_bg), mean(chip_tmp1_bg), col='black', pch=16)
abline(0,1,col='red')
dev.off()


### plot qtPCR vs chipseq normed (whole genome)
for (i in c(1:4)){
print(i)
chip_tmp1_pk = chipseq_od_mat_all_pk[sample_id,2] < 0.05
chip_tmp2_pk = chipseq_od_mat_all_pk[sample_id,i] < 0.05
cpk = chip_tmp1_pk*chip_tmp2_pk > 0
cbg = chip_tmp1_pk+chip_tmp2_pk == 0

A1 = test_AB_mat[1,2]
B1 = test_AB_mat[2,2]
A2 = test_AB_mat[1,i]
B2 = test_AB_mat[2,i]
print('A2-B2:')
print(A2)
print(B2)
png(paste('test.norm.chipseq_qtpcr.s3norm.compare0HR.', toString(i), '.png', sep=''))
chip_tmp1 = log2(A1 * (chipseq_od_mat_all[sample_id,2]+small_num)^B1)
chip_tmp2 =  log2(A2 * (chipseq_od_mat_all[sample_id,i]+small_num)^B2)
chip_tmp1_pk = chipseq_od_mat_all_pk[sample_id,2] < 0.05
chip_tmp2_pk = chipseq_od_mat_all_pk[sample_id,i] < 0.05
chip_tmp1_qtpcr_region = log2(A1 * (chip_fihserp[,2]+small_num)^B1)
chip_tmp2_qtpcr_region =  B2 * log2((chip_fihserp[,i]+small_num)) + A2
chip_tmp2_qtpcr_region_od =  log2((chip_fihserp[,i]+small_num))
cpk = chip_tmp1_pk*chip_tmp2_pk > 0
cbg = chip_tmp1_pk+chip_tmp2_pk == 0
log2qtpcr_norm_tmp1 = log2qtpcr_norm[,2]
log2qtpcr_norm_tmp2 = log2qtpcr_norm[,i]
print('correlation')
print(cor(log2qtpcr_norm_tmp2, chip_tmp2_qtpcr_region))
print(cor(log2qtpcr_norm_tmp2, chip_tmp2_qtpcr_region_od))
print(cor(chip_tmp2_qtpcr_region, chip_tmp2_qtpcr_region_od))
print((chip_tmp2_qtpcr_region))
print((chip_tmp2_qtpcr_region_od))
stable_pr2 = cor(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic'])
plot(chip_tmp1, chip_tmp2, xlim=c(-3, 9), ylim=c(-3, 9), col='gray', main=toString(round(stable_pr2, 3)))
#points(chip_tmp1[cpk], chip_tmp2[cpk], col='orangered')
#plot(log2qtpcr_norm_tmp1, log2qtpcr_norm_tmp2, xlim=c(-3, 9), ylim=c(-3, 9))
points(chip_tmp1_qtpcr_region[regions_name=='Strong-Dynamic'], chip_tmp2_qtpcr_region[regions_name=='Strong-Dynamic'], col='blue')
points(chip_tmp1_qtpcr_region[regions_name=='Weak'], chip_tmp2_qtpcr_region[regions_name=='Weak'], col='green')
points(chip_tmp1_qtpcr_region[regions_name=='Strong-Stable'], chip_tmp2_qtpcr_region[regions_name=='Strong-Stable'], col='red')
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name=='Strong-Dynamic'], col='blue', pch=16)
points(log2qtpcr_norm_tmp1[regions_name=='Weak'], log2qtpcr_norm_tmp2[regions_name=='Weak'], col='green', pch=16)
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable'], log2qtpcr_norm_tmp2[regions_name=='Strong-Stable'], col='red', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable']), mean(log2qtpcr_norm_tmp2[regions_name=='Strong-Stable']), col='red', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Weak']), mean(log2qtpcr_norm_tmp2[regions_name=='Weak']), col='green', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic']), mean(log2qtpcr_norm_tmp2[regions_name=='Strong-Dynamic']), col='blue', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name!='Strong-Dynamic']), mean(log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']), col='deepskyblue', pch=16)
points(mean(chip_tmp1), mean(chip_tmp2), col='deepskyblue')
points(mean(chip_tmp1[cpk]), mean(chip_tmp2[cpk]), col='orangered')
points(mean(chip_tmp1[cbg]), mean(chip_tmp2[cbg]), col='black')
abline(0,1,col='black')
dev.off()
}



















### plot qtPCR vs chipseq (qtPCR regions)
for (i in c(1:4)){
png(paste('test.norm.chipseq_qtpcr.qtpcr_region.', toString(i), '.png', sep=''))
chip_tmp1 = log2chipseq[,i]
chip_tmp1_all = log2(chipseq_od_mat_all[sample_id,i]+small_num)
log2qtpcr_norm_tmp1 = log2qtpcr_norm[,i]
plot(log2qtpcr_norm_tmp1, chip_tmp1, xlim=c(-3, 9), ylim=c(-3, 9), col='gray')
#plot(log2qtpcr_norm_tmp1, log2qtpcr_norm_tmp2, xlim=c(-3, 9), ylim=c(-3, 9))
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic'], chip_tmp1[regions_name=='Strong-Dynamic'], col='blue')
points(log2qtpcr_norm_tmp1[regions_name=='Weak'], chip_tmp1[regions_name=='Weak'], col='green')
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable'], chip_tmp1[regions_name=='Strong-Stable'], col='red')
points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable']), mean(chip_tmp1[regions_name=='Strong-Stable']), col='red', pch=16)
points(mean(log2qtpcr_norm_tmp1[regions_name=='Weak']), mean(chip_tmp1[regions_name=='Weak']), col='green', pch=16)
points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic']), mean(chip_tmp1[regions_name=='Strong-Dynamic']), col='blue', pch=16)
points(mean(log2qtpcr_norm_tmp1[regions_name!='Strong-Dynamic']), mean(chip_tmp1[regions_name!='Strong-Dynamic']), col='deepskyblue', pch=16)
points(mean(log2qtpcr_norm_tmp1), mean(chip_tmp1), col='black', pch=16)
points(mean(chip_tmp1_all), mean(chip_tmp1_all), col='black', pch=16)
abline(0,1,col='red')
dev.off()
}





### plot qtPCR vs chipseq (whole genome)
for (i in c(1:4)){
print(i)
A1 = test_AB_mat[1,2]
B1 = test_AB_mat[2,2]
A2 = test_AB_mat[1,i]
B2 = test_AB_mat[2,i]
png(paste('test.norm.chipseq_qtpcr.s3norm.compare0HR.', toString(i), '.png', sep=''))
chip_tmp1 = B1 * log2(chipseq_od_mat_all[sample_id,2]+small_num) + A1
chip_tmp2 = B2 * log2(chipseq_od_mat_all[sample_id,i]+small_num) + A2
chip_tmp1_qtpcr_region = B1 * log2chipseq[,2] + A1
chip_tmp2_qtpcr_region = B2 * log2chipseq[,i] + A2
print(A2)
print(B2)

chip_tmp1_pk = chipseq_od_mat_all_pk[sample_id,2] < 0.05
chip_tmp2_pk = chipseq_od_mat_all_pk[sample_id,i] < 0.05
cpk = chip_tmp1_pk*chip_tmp2_pk > 0
cbg = chip_tmp1_pk+chip_tmp2_pk == 0
log2qtpcr_norm_tmp1 = log2qtpcr_norm[,2]
log2qtpcr_norm_tmp2 = log2qtpcr_norm[,i]
print(cor.test(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']))
cor_test = cor.test(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic'])
stable_pr2 = cor(2^chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'], 2^log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic'])
#stable_pr2 = cor(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic'])
stable_pr2 = 1 - sum((log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']-chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic'])^2) / sum((log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']-mean(log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']))^2) 

plot(chip_tmp1, chip_tmp2, xlim=c(-3, 9), ylim=c(-3, 9), col='gray', main=paste(toString(round(cor_test$estimate, 3)), toString(cor_test$p.value), sep=' ') )
#points(chip_tmp1[cpk], chip_tmp2[cpk], col='orangered')
#plot(log2qtpcr_norm_tmp1, log2qtpcr_norm_tmp2, xlim=c(-3, 9), ylim=c(-3, 9))
### plot chip-seq signal at qtPCR regions
points(chip_tmp1_qtpcr_region[regions_name=='Strong-Dynamic'], chip_tmp2_qtpcr_region[regions_name=='Strong-Dynamic'], col='blue')
points(chip_tmp1_qtpcr_region[regions_name=='Weak'], chip_tmp2_qtpcr_region[regions_name=='Weak'], col='green')
points(chip_tmp1_qtpcr_region[regions_name=='Strong-Stable'], chip_tmp2_qtpcr_region[regions_name=='Strong-Stable'], col='red')
lm_stable_chip = lm(chip_tmp2_qtpcr_region[regions_name!='Strong-Dynamic']~chip_tmp1_qtpcr_region[regions_name!='Strong-Dynamic'])
abline(lm_stable_chip, col="blue")
### plot qtPCR signal at qtPCR regions
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic'], log2qtpcr_norm_tmp2[regions_name=='Strong-Dynamic'], col='blue', pch=16)
points(log2qtpcr_norm_tmp1[regions_name=='Weak'], log2qtpcr_norm_tmp2[regions_name=='Weak'], col='green', pch=16)
points(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable'], log2qtpcr_norm_tmp2[regions_name=='Strong-Stable'], col='red', pch=16)
lm_stable_qtpcr = lm(log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']~log2qtpcr_norm_tmp1[regions_name!='Strong-Dynamic'])
abline(lm_stable_qtpcr, col="red")
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Stable']), mean(log2qtpcr_norm_tmp2[regions_name=='Strong-Stable']), col='red', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Weak']), mean(log2qtpcr_norm_tmp2[regions_name=='Weak']), col='green', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name=='Strong-Dynamic']), mean(log2qtpcr_norm_tmp2[regions_name=='Strong-Dynamic']), col='blue', pch=16)
#points(mean(log2qtpcr_norm_tmp1[regions_name!='Strong-Dynamic']), mean(log2qtpcr_norm_tmp2[regions_name!='Strong-Dynamic']), col='deepskyblue', pch=16)
points(mean(chip_tmp1), mean(chip_tmp2), col='deepskyblue')
points(mean(chip_tmp1[cpk]), mean(chip_tmp2[cpk]), col='orangered')
points(mean(chip_tmp1[cbg]), mean(chip_tmp2[cbg]), col='black')
abline(0,1,col='black')
dev.off()
}


