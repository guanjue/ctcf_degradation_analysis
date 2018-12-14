library(MASS)

#######################################
### FOR 5T7D
#######################################
print('FOR 5T7D')
all_sample_input_table = read.table('bam_bw_bedgraph/file_list.all.txt', header=F)[c(1,2,16:19, 3,4, 20:24),]
IP_list = all_sample_input_table[,1]
CTRL_list = all_sample_input_table[,2]
TP_list = apply(all_sample_input_table, 1, function(x) unlist(strsplit(x[1], split='_'))[2])
ID_list = apply(all_sample_input_table, 1, function(x) unlist(strsplit(x[1], split='_'))[1])

### initialize matrices
qPCR_expand_mat = c()
d_IP_raw_mat = c()
d_IP_ts_mat = c()
d_IP_s3_mat = c()
sf_mat = c()
sf_lm_mat = c()
sf_lm_mat_A = c()
cor_vec_ts = c()
cor_vec_s3 = c()
name_mat = c()

method_list = c('raw', 'sub', 'log2fc', 'fc0', 'fc')
method_name = c('RAW', 'Sub', 'log2FC', 'FC', 'FC_MACS', 'log2_FC_MACS')
method_name_s3 = c('RAW_LM', 'Sub_LM', 'log2FC_LM', 'FC_LM', 'FC_MACS_LM', 'log2_FC_MACS_LM')

qPCR_len = read.table('qtPCR_sig.5T7D.len.txt', header=F)
used_id_len = qPCR_len>=400
qPCR_pk_id = read.table('qtPCR_sig.5T7D.ids.txt', header=F)[used_id_len,]

###########
s3norm_no0 = function(ref, tar){
	used_id = ref!=0
	fit = lm(ref[used_id]~tar[used_id])
	used_AB = fit$coefficients
	return(used_AB)
}
###########
###########
s3norm = function(ref, tar){
	fit = lm(ref~tar-1)
	used_AB = fit$coefficients
	return(used_AB)
}
###########
###########
s3normlog2 = function(ref, tar){
	ref_log2 = log2(ref+1)
	tar_log2 = log2(tar+1)
	fit = lm(ref_log2~tar_log2)
	used_AB = fit$coefficients
	return(used_AB)
}
###########
R2 = function(d_IP_ts, qPCR_expand, used_id){
	r2 = 1-sum((qPCR_expand-d_IP_ts)[used_id]^2) / sum((qPCR_expand-mean(qPCR_expand))[used_id]^2)
	return(r2)
}
###########

###########
R2_all = function(d_IP_ts, qPCR_expand){
        r2 = 1-sum((qPCR_expand-d_IP_ts)^2) / sum((qPCR_expand-mean(qPCR_expand))^2)
        return(r2)
}
###########

qPCR_len = read.table('qtPCR_sig.5T7D.len.txt', header=F)

### get qPCR norm results
for (i in c(1:length(IP_list))){
	print(IP_list[i])
	qPCR = scan(paste('qtPCR_sig.5T7D.', toString(TP_list[i]), '.txt', sep=''))[used_id_len]
	used_id = qPCR!=-10
	qPCR_expand = qPCR * 100
	for (j in c(1:length(method_list))){
		d_IP = scan(paste('qtPCR_regions_signal_tab/', toString(IP_list[i]), '.qtPCR_regions.sorted.signal.5T7D.', method_list[j], '.tab', sep=''))[used_id_len]
		### get TSnorm scale factor
		sf_tmp = 1/mean(d_IP[used_id]) * mean(qPCR_expand[used_id])
		sf_s3_tmp = s3norm((qPCR_expand[used_id]), (d_IP[used_id]))
		### get TSnorm 
		d_IP_ts = d_IP * sf_tmp
		d_IP_s3 = d_IP * sf_s3_tmp
		d_IP_s3[d_IP_s3<0]=0
		### get correlation
		cor_tmp_ts = 1-sum((qPCR_expand-d_IP_ts)[used_id]^2) / sum((qPCR_expand-mean(qPCR_expand))[used_id]^2) #cor(qPCR_expand, d_IP, method='pearson')
		cor_tmp_s3 = 1-sum((qPCR_expand-d_IP_s3)[used_id]^2) / sum((qPCR_expand-mean(qPCR_expand))[used_id]^2) #cor(qPCR_expand, d_IP_ts, method='pearson')
		### get names
		name_tmp = paste(ID_list[i], method_list[j], TP_list[i], sep='_')
		###### collect results
		qPCR_expand_mat = cbind(qPCR_expand_mat, qPCR_expand)
		d_IP_raw_mat = cbind(d_IP_raw_mat, d_IP)
		d_IP_ts_mat = cbind(d_IP_ts_mat, d_IP_ts)
		d_IP_s3_mat = cbind(d_IP_s3_mat, d_IP_s3)
		sf_mat = cbind(sf_mat, sf_tmp)
		sf_lm_mat = cbind(sf_lm_mat, sf_s3_tmp)
		cor_vec_ts = cbind(cor_vec_ts, cor_tmp_ts)
		cor_vec_s3 = cbind(cor_vec_s3, cor_tmp_s3)
		name_mat = cbind(name_mat, name_tmp)
	}
	### get log2FC_MACS
	d_IP0 = scan(paste('qtPCR_regions_signal_tab/', toString(IP_list[i]), '.qtPCR_regions.sorted.signal.5T7D.', 'fc', '.tab', sep=''))[used_id_len]
	d_IP = log2(d_IP0+1)
	### get TSnorm scale factor
	sf_tmp = 1/mean(d_IP[used_id]) * mean(qPCR_expand[used_id])
	sf_s3_tmp = s3norm((qPCR_expand[used_id]), (d_IP[used_id]))
	### get TSnorm 
	d_IP_ts = d_IP * sf_tmp
	d_IP_s3 = d_IP * sf_s3_tmp
	d_IP_s3[d_IP_s3<0]=0
	#d_IP_s3 = 2^(log2(d_IP) * sf_s3_tmp[2] + sf_s3_tmp[1])
	### get correlation
	cor_tmp_ts = 1-sum((qPCR_expand-d_IP_ts)[used_id]^2) / sum((qPCR_expand-mean(qPCR_expand))[used_id]^2) #cor(qPCR_expand, d_IP, method='pearson')
	cor_tmp_s3 = 1-sum((qPCR_expand-d_IP_s3)[used_id]^2) / sum((qPCR_expand-mean(qPCR_expand))[used_id]^2) #cor(qPCR_expand, d_IP_ts, method='pearson')
	### get names
	name_tmp = paste(ID_list[i], 'log2_FC_MACS', TP_list[i], sep='_')
	###### collect results
	qPCR_expand_mat = cbind(qPCR_expand_mat, qPCR_expand)
	d_IP_raw_mat = cbind(d_IP_raw_mat, d_IP)
	d_IP_ts_mat = cbind(d_IP_ts_mat, d_IP_ts)
	d_IP_s3_mat = cbind(d_IP_s3_mat, d_IP_s3)
	sf_mat = cbind(sf_mat, sf_tmp)
	sf_lm_mat = cbind(sf_lm_mat, sf_s3_tmp)
	cor_vec_ts = cbind(cor_vec_ts, cor_tmp_ts)
	cor_vec_s3 = cbind(cor_vec_s3, cor_tmp_s3)
	name_mat = cbind(name_mat, name_tmp)
}

print('get colnames!!!')
### get colnames
colnames(qPCR_expand_mat) = name_mat
colnames(d_IP_raw_mat) = name_mat
colnames(d_IP_ts_mat) = name_mat
colnames(d_IP_s3_mat) = name_mat
colnames(sf_mat) = name_mat
colnames(sf_lm_mat) = name_mat
colnames(cor_vec_ts) = name_mat
colnames(cor_vec_s3) = name_mat
print('get rownames!!!')
### get rownames
print('get rownames!!!')
print(dim(qPCR_expand_mat))
print(length(qPCR_pk_id))
rownames(qPCR_expand_mat) = qPCR_pk_id
print('get rownames!!!')
rownames(d_IP_raw_mat) = qPCR_pk_id
rownames(d_IP_ts_mat) = qPCR_pk_id
rownames(d_IP_s3_mat) = qPCR_pk_id
rownames(sf_mat) = 'SF'
rownames(sf_lm_mat) = c('B')
rownames(cor_vec_ts) = 'COR'
rownames(cor_vec_s3) = 'COR'
print('reshape sf_mat & sf_mat!!!')
### reshape sf_mat & sf_mat
sf_mat = matrix(sf_mat, nrow = length(method_name), byrow = FALSE)
colnames(sf_mat) = ID_list
rownames(sf_mat) = method_name

sf_lm_mat = matrix(sf_lm_mat, nrow = length(method_name), byrow = FALSE)
colnames(sf_lm_mat) = ID_list
rownames(sf_lm_mat) = method_name

cor_mat_ts = matrix(cor_vec_ts, nrow = length(method_name), byrow = FALSE)
colnames(cor_mat_ts) = ID_list
rownames(cor_mat_ts) = method_name

cor_mat_s3 = matrix(cor_vec_s3, nrow = length(method_name), byrow = FALSE)
colnames(cor_mat_s3) = ID_list
rownames(cor_mat_s3) = method_name


cor_mat_ts_s3 = rbind(cor_mat_ts, cor_mat_s3)
### output
dir.create('qPCR_cor_5T7D')

### plot correlation boxplot
pdf('qPCR_cor_5T7D/correlation_boxplot_ts_s3.pdf')
cor_mat_t0 = t(cor_mat_ts_s3)
cor_mat_t = c()
method_2 = c()
for (i in c(1:(dim(cor_mat_t0)[2]/2))){
	cor_mat_t = cbind(cor_mat_t, cor_mat_t0[,i])
	cor_mat_t = cbind(cor_mat_t, cor_mat_t0[,i+(dim(cor_mat_t0)[2]/2)])
	method_2 = c(method_2, method_name[i], method_name_s3[i])
}
boxplot(cor_mat_t, ylim=c(-1, 1) ,outline=FALSE, xaxt="n", xlab='Different Methods', main='R2 = 1 - SSE/SSTO', ylab='R2')
axis(1, at=c(1:dim(cor_mat_t)[2]), labels=FALSE)
text(cex=1, x=c(1:dim(cor_mat_t)[2]), y=-0.5, method_2, xpd=TRUE, srt=90)
for (i in c(1:dim(cor_mat_t)[1])){
	lines(c(1:(dim(cor_mat_t)[2])), cor_mat_t[i,], col=rgb(255/255,0/255,0/255,alpha=0.5) )
	points(c(1:(dim(cor_mat_t)[2])), cor_mat_t[i,], col='black', pch=20 )
}
dev.off()

### plot scatterplot
#plot_lim = max(abs(cbind(d_IP_ts_mat, qPCR_expand_mat)))
plot_lim = 35
for (i in c(1: dim(d_IP_ts_mat)[2])){
	output_fig_name = paste('qPCR_cor_5T7D/', colnames(d_IP_ts_mat)[i], '.scatterplot.pdf', sep='')
	pdf(output_fig_name, width=9, height=5)
	par(mfrow=c(1,2))
	#plot(qPCR_expand_mat[,i], d_IP_raw_mat[,i], pch = 20, col='black', xlim=c(0, plot_lim), ylim=c(0, plot_lim), xlab='qPCR*100', ylab='TSnorm-ChIP-seq', main = paste(colnames(d_IP_ts_mat)[i], ': ', toString(round(cor_vec_ts[i], 2)), sep='') )
	#abline(0,1, col='red', lwd = 1.5)
	plot(qPCR_expand_mat[,i], d_IP_ts_mat[,i], pch = 20, col='black', xlim=c(0, plot_lim), ylim=c(0, plot_lim), xlab='qPCR*100', ylab='TSnorm-ChIP-seq', main = paste(colnames(d_IP_ts_mat)[i], ': ', toString(round(cor_vec_ts[i], 2)), sep='') )
	abline(0,1, col='red', lwd = 1.5)
	plot(qPCR_expand_mat[,i], d_IP_s3_mat[,i], pch = 20, col='black', xlim=c(0, plot_lim), ylim=c(0, plot_lim), xlab='qPCR*100', ylab='TSnorm-ChIP-seq', main = paste(colnames(d_IP_ts_mat)[i], ': ', toString(round(cor_vec_s3[i], 2)), sep='') )
	abline(0,1, col='red', lwd = 1.5)
	dev.off()
}

### plot barplot
get_barplot_fig = function(qPCR_expand_mat, method_name, k){
	### get qPCR
	qPCR_expand_mat_tmp = c()
	id_tmp = c()
	for (i in c(1:(dim(qPCR_expand_mat)[2]/length(method_name))) ){
		qPCR_expand_vec_tmp = qPCR_expand_mat[,k+(i-1)*length(method_name)]
		qPCR_expand_mat_tmp = cbind(qPCR_expand_mat_tmp, qPCR_expand_vec_tmp)
		id_tmp[i] = paste(unlist(strsplit(colnames(qPCR_expand_mat)[k+(i-1)*length(method_name)], '_'))[1], unlist(strsplit(colnames(qPCR_expand_mat)[k+(i-1)*length(method_name)], '_'))[3], sep='_')
	}
	### get colnames
	colnames(qPCR_expand_mat_tmp) = id_tmp
	###
	qPCR_expand_mat_TP = c()
	for (TPtmp in unique(TP_list)){
	qPCR_expand_mat_TP_tmp = apply(qPCR_expand_mat_tmp[,TP_list==TPtmp], 1, mean)
	qPCR_expand_mat_TP = cbind(qPCR_expand_mat_TP, qPCR_expand_mat_TP_tmp)
	}
	colnames(qPCR_expand_mat_TP) = unique(TP_list)
	return(qPCR_expand_mat_TP)
}

###
colours = c("red","orange",'yellow','green','blue','purple')

pdf('qPCR_cor_5T7D/6tp_bar.qPCR_5T7D.pdf', width=12, height=12)
par(mfrow=c(4,1))

qPCR_expand_mat_TP = get_barplot_fig(qPCR_expand_mat, method_name, 1)
barplot(t(qPCR_expand_mat_TP), beside = TRUE, col=colours, xaxt="n", main = 'qPCR')
axis(1, at=c(1:dim(t(qPCR_expand_mat_TP))[2])*7, labels=FALSE)
text(cex=1, x=3.5+(c(1:dim(t(qPCR_expand_mat_TP))[2])-1)*7, y=-5, rownames(qPCR_expand_mat_TP), xpd=TRUE, srt=45)
box()

d_IP_raw_mat_TP = get_barplot_fig(d_IP_raw_mat, method_name, 5)
main_name = 'RAW'
for (i in c(1:dim(d_IP_raw_mat_TP)[2])){
cor_1 = R2(d_IP_raw_mat_TP[,i], qPCR_expand_mat_TP[,i], qPCR_expand_mat_TP[,i]!=0)
main_name = paste(main_name, toString(round(cor_1,2)), sep=':')
}
barplot(t(d_IP_raw_mat_TP), beside = TRUE, col=colours, xaxt="n", main = main_name)
axis(1, at=c(1:dim(t(d_IP_raw_mat_TP))[2])*7, labels=FALSE)
text(cex=1, x=3.5+(c(1:dim(t(d_IP_raw_mat_TP))[2])-1)*7, y=-40, rownames(d_IP_raw_mat_TP), xpd=TRUE, srt=45)
box()

d_IP_ts_mat_TP = get_barplot_fig(d_IP_ts_mat, method_name, 5)
main_name = 'TSnorm'
for (i in c(1:dim(d_IP_ts_mat_TP)[2])){
cor_1 = R2(d_IP_ts_mat_TP[,i], qPCR_expand_mat_TP[,i], qPCR_expand_mat_TP[,i]!=0)
main_name = paste(main_name, toString(round(cor_1,2)), sep=':')
}
barplot(t(d_IP_ts_mat_TP), beside = TRUE, col=colours, xaxt="n", main = main_name)
axis(1, at=c(1:dim(t(d_IP_ts_mat_TP))[2])*7, labels=FALSE)
text(cex=1, x=3.5+(c(1:dim(t(d_IP_ts_mat_TP))[2])-1)*7, y=-5, rownames(d_IP_ts_mat_TP), xpd=TRUE, srt=45)
box()

d_IP_s3_mat_TP = get_barplot_fig(d_IP_s3_mat, method_name, 5)
print('!!!!!')
print(dim(d_IP_s3_mat_TP))
main_name = 'LM'
for (i in c(1:dim(d_IP_s3_mat_TP)[2])){
cor_1 = R2(d_IP_s3_mat_TP[,i], qPCR_expand_mat_TP[,i], qPCR_expand_mat_TP[,i]!=0)
main_name = paste(main_name, toString(round(cor_1,2)), sep=':')
}
barplot(t(d_IP_s3_mat_TP), beside = TRUE, col=colours, xaxt="n", main = main_name)
axis(1, at=c(1:dim(t(d_IP_s3_mat_TP))[2])*7, labels=FALSE)
text(cex=1, x=3.5+(c(1:dim(t(d_IP_s3_mat_TP))[2])-1)*7, y=-7, rownames(d_IP_s3_mat_TP), xpd=TRUE, srt=45)
box()

dev.off()

### write all scatterplot
pdf('qPCR_cor_5T7D/d_IP_s3_mat.pdf', width=12, height=7)
par(mfrow=c(1,2))
d_IP_ts_mat_TP_vec = as.vector(d_IP_ts_mat_TP)
d_IP_s3_mat_TP_vec = as.vector(d_IP_s3_mat_TP)
qPCR_expand_mat_TP_vec = as.vector(qPCR_expand_mat_TP)

R_ts = cor(d_IP_ts_mat_TP_vec[qPCR_expand_mat_TP_vec!=0], qPCR_expand_mat_TP_vec[qPCR_expand_mat_TP_vec!=0])
R2_ts = R2_all(d_IP_ts_mat_TP_vec[qPCR_expand_mat_TP_vec!=0], qPCR_expand_mat_TP_vec[qPCR_expand_mat_TP_vec!=0])

R_s3 = cor(d_IP_s3_mat_TP_vec[qPCR_expand_mat_TP_vec!=0], qPCR_expand_mat_TP_vec[qPCR_expand_mat_TP_vec!=0])
R2_s3 = R2_all(d_IP_s3_mat_TP_vec[qPCR_expand_mat_TP_vec!=0], qPCR_expand_mat_TP_vec[qPCR_expand_mat_TP_vec!=0])

main_name_ts = paste('R: ', round(R_ts,2), '; ', 'R2: ', round(R2_ts,2),  sep='')
main_name_s3 = paste('R: ', round(R_s3,2), '; ', 'R2: ', round(R2_s3,2),  sep='')
plot(qPCR_expand_mat_TP_vec, d_IP_ts_mat_TP_vec, xlim=c(0, plot_lim), ylim=c(0, plot_lim), main=main_name_ts, xlab = 'qPCR*100', ylab = 'TSnorm')
abline(0,1,col='red')
plot(qPCR_expand_mat_TP_vec, d_IP_s3_mat_TP_vec, xlim=c(0, plot_lim), ylim=c(0, plot_lim), main=main_name_s3, xlab = 'qPCR*100', ylab = 'LM')
abline(0,1,col='red')
dev.off()

### write scale factors
write.table(sf_mat, 'TSnorm_SF.5T7D.txt', col.names=TRUE, row.names=TRUE, quote=FALSE, sep='\t')
write.table(sf_lm_mat, 's3norm_SF.5T7D.txt', col.names=TRUE, row.names=TRUE, quote=FALSE, sep='\t')


