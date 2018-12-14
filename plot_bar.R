cor_mat = c()
qPCR_mat = c()
d0_mat = c()
d0_ts_mat = c()


qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.6TP.0A.txt')
d1 = scan('qtPCR_regions_signal_tab/1579_0A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d2 = scan('qtPCR_regions_signal_tab/1580_0A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d0 = (d1+d2)/2
qPCR_mat = cbind(qPCR_mat, qPCR)
d0_mat = cbind(d0_mat, d0)
d0_ts_mat = cbind(d0_ts_mat, d0/mean(d0)*mean(qPCR))
cor_mat = cbind(cor_mat, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.6TP.4A.txt')
d1 = scan('qtPCR_regions_signal_tab/1582_4A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d2 = scan('qtPCR_regions_signal_tab/1583_4A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d0 = (d1+d2)/2
qPCR_mat = cbind(qPCR_mat, qPCR)
d0_mat = cbind(d0_mat, d0)
d0_ts_mat = cbind(d0_ts_mat, d0/mean(d0)*mean(qPCR))
cor_mat = cbind(cor_mat, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.6TP.6A.txt')
d1 = scan('qtPCR_regions_signal_tab/1585_6A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d2 = scan('qtPCR_regions_signal_tab/1586_6A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d0 = (d1+d2)/2
qPCR_mat = cbind(qPCR_mat, qPCR)
d0_mat = cbind(d0_mat, d0)
d0_ts_mat = cbind(d0_ts_mat, d0/mean(d0)*mean(qPCR))
cor_mat = cbind(cor_mat, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.6TP.12A.txt')
d1 = scan('qtPCR_regions_signal_tab/1849_12A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d2 = scan('qtPCR_regions_signal_tab/1850_12A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d0 = (d1+d2)/2
qPCR_mat = cbind(qPCR_mat, qPCR)
d0_mat = cbind(d0_mat, d0)
d0_ts_mat = cbind(d0_ts_mat, d0/mean(d0)*mean(qPCR))
cor_mat = cbind(cor_mat, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.6TP.18A.txt')
d1 = scan('qtPCR_regions_signal_tab/1861_18A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d2 = scan('qtPCR_regions_signal_tab/1886_18A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d0 = (d1+d2)/2
qPCR_mat = cbind(qPCR_mat, qPCR)
d0_mat = cbind(d0_mat, d0)
d0_ts_mat = cbind(d0_ts_mat, d0/mean(d0)*mean(qPCR))
cor_mat = cbind(cor_mat, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.6TP.24A.txt')
d1 = scan('qtPCR_regions_signal_tab/1863_24A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d2 = scan('qtPCR_regions_signal_tab/1880_24A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.6TP.tab')
d0 = (d1+d2)/2
qPCR_mat = cbind(qPCR_mat, qPCR)
d0_mat = cbind(d0_mat, d0)
d0_ts_mat = cbind(d0_ts_mat, d0/mean(d0)*mean(qPCR))
cor_mat = cbind(cor_mat, cor(qPCR, d1+d2, method='pearson'))



colours = c("red","orange",'yellow','green','blue', 'purple')

pdf('6tp_bar.qPCR.pdf', width=14, height=21)
par(mfrow=c(3,1))
barplot(t(qPCR_mat),beside = TRUE, col=colours)
barplot(t(d0_mat),beside = TRUE, col=colours)
barplot(t(d0_ts_mat),beside = TRUE, col=colours)
dev.off()










cor_mat_5T7D = c()
qPCR_mat_5T7D = c()
d0_mat_5T7D = c()
d0_ts_mat_5T7D = c()

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.5T7D.0A.txt')
d1 = scan('qtPCR_regions_signal_tab/1579_0A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d2 = scan('qtPCR_regions_signal_tab/1580_0A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d0 = (d1+d2)/2
qPCR_mat_5T7D = cbind(qPCR_mat_5T7D, qPCR)
d0_mat_5T7D = cbind(d0_mat_5T7D, d0)
d0_ts_mat_5T7D = cbind(d0_ts_mat_5T7D, d0/mean(d0)*mean(qPCR))
cor_mat_5T7D = cbind(cor_mat_5T7D, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.5T7D.5T.txt')
d1 = scan('qtPCR_regions_signal_tab/1917_5T_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d2 = scan('qtPCR_regions_signal_tab/1918_5T_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d0 = (d1+d2)/2
qPCR_mat_5T7D = cbind(qPCR_mat_5T7D, qPCR)
d0_mat_5T7D = cbind(d0_mat_5T7D, d0)
d0_ts_mat_5T7D = cbind(d0_ts_mat_5T7D, d0/mean(d0)*mean(qPCR))
cor_mat_5T7D = cbind(cor_mat_5T7D, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.5T7D.7D.txt')
d1 = scan('qtPCR_regions_signal_tab/1920_7D_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d2 = scan('qtPCR_regions_signal_tab/1921_7D_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d0 = (d1+d2)/2
qPCR_mat_5T7D = cbind(qPCR_mat_5T7D, qPCR)
d0_mat_5T7D = cbind(d0_mat_5T7D, d0)
d0_ts_mat_5T7D = cbind(d0_ts_mat_5T7D, d0/mean(d0)*mean(qPCR))
cor_mat_5T7D = cbind(cor_mat_5T7D, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.5T7D.4A.txt')
d1 = scan('qtPCR_regions_signal_tab/1582_4A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d2 = scan('qtPCR_regions_signal_tab/1583_4A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d0 = (d1+d2)/2
qPCR_mat_5T7D = cbind(qPCR_mat_5T7D, qPCR)
d0_mat_5T7D = cbind(d0_mat_5T7D, d0)
d0_ts_mat_5T7D = cbind(d0_ts_mat_5T7D, d0/mean(d0)*mean(qPCR))
cor_mat_5T7D = cbind(cor_mat_5T7D, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.5T7D.5T4A.txt')
d1 = scan('qtPCR_regions_signal_tab/1859_5T4A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d2 = scan('qtPCR_regions_signal_tab/1881_5T4A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d0 = (d1+d2)/2
qPCR_mat_5T7D = cbind(qPCR_mat_5T7D, qPCR)
d0_mat_5T7D = cbind(d0_mat_5T7D, d0)
d0_ts_mat_5T7D = cbind(d0_ts_mat_5T7D, d0/mean(d0)*mean(qPCR))
cor_mat_5T7D = cbind(cor_mat_5T7D, cor(qPCR, d1+d2, method='pearson'))

qPCR = scan('/storage/home/gzx103/group/projects/ctcf_auxin/qtPCR_sig.5T7D.7D4A.txt')
d1 = scan('qtPCR_regions_signal_tab/1857_7D4A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d2 = scan('qtPCR_regions_signal_tab/1882_7D4A_CTCF_mm9_duprm.bam.qtPCR_regions.sorted.signal.5T7D.tab')
d0 = (d1+d2)/2
qPCR_mat_5T7D = cbind(qPCR_mat_5T7D, qPCR)
d0_mat_5T7D = cbind(d0_mat_5T7D, d0)
d0_ts_mat_5T7D = cbind(d0_ts_mat_5T7D, d0/mean(d0)*mean(qPCR))
cor_mat_5T7D = cbind(cor_mat_5T7D, cor(qPCR, d1+d2, method='pearson'))


pdf('5t7d_bar.qPCR.pdf', width=14, height=21)
par(mfrow=c(3,1))
barplot(t(qPCR_mat_5T7D),beside = TRUE, col=colours)
barplot(t(d0_mat_5T7D),beside = TRUE, col=colours)
barplot(t(d0_ts_mat_5T7D),beside = TRUE, col=colours)
dev.off()




