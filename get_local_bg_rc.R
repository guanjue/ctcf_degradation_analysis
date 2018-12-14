args = commandArgs(trailingOnly=TRUE)
wg_tab = args[1]
input_rc_1kb = args[2]
input_rc_5kb = args[3]
input_rc_10kb = args[4]
output = args[5]

wg_tab_mat = read.table(wg_tab)
wg_mean_input_sig = sum(wg_tab_mat[,4]) / sum(wg_tab_mat[,2])

input_rc_1kb_sig = scan(input_rc_1kb)
input_rc_5kb_sig = scan(input_rc_5kb)
input_rc_10kb_sig = scan(input_rc_10kb)

wg_mean_input_rc_sig = rep(wg_mean_input_sig, length(input_rc_10kb_sig))

input_sig_mat = cbind(wg_mean_input_rc_sig, input_rc_1kb_sig, input_rc_5kb_sig, input_rc_10kb_sig)

input_sig_mat_max = apply(input_sig_mat, 1, max)

print(dim(input_sig_mat_max))

write.table(input_sig_mat_max, output, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')



