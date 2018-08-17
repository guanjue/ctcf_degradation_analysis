#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc
module load python/2.7.14-anaconda5.0.1

cd /storage/home/gzx103/scratch/ctcf

input_dir=/storage/home/gzx103/scratch/ctcf/
script_dir=/storage/home/gzx103/group/software/PKnorm/pknorm_scripts/
working_dir=/storage/home/gzx103/scratch/ctcf/

while read files
do
	ct=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	ctrl=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}' | awk -F '_' -v OFS='\t' '{print $1}')
	rep1=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $3}' | awk -F '_' -v OFS='\t' '{print $1}')
	rep2=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $4}' | awk -F '_' -v OFS='\t' '{print $1}')
	echo $ct
	echo $ctrl
	echo $rep1
	echo $rep2
	cut -f4 $ctrl'_200bp_rc.txt' > $ctrl'_rc_all.txt'
	cut -f4 $rep1'_200bp_rc.txt' > $rep1'_rc_all.txt'
	cut -f4 $rep2'_200bp_rc.txt' > $rep2'_rc_all.txt'
	cut -f4 $ctrl'_qtpcr_rc.txt' >> $ctrl'_rc_all.txt'
	cut -f4 $rep1'_qtpcr_rc.txt' >> $rep1'_rc_all.txt'
	cut -f4 $rep2'_qtpcr_rc.txt' >> $rep2'_rc_all.txt'
	time Rscript $script_dir'negative_binomial_p_2r_bgadj_bayes.R' $rep1'_rc_all.txt' $input_dir $ctrl'_rc_all.txt' $input_dir $ct'.rep1'
	time Rscript $script_dir'negative_binomial_p_2r_bgadj_bayes.R' $rep2'_rc_all.txt' $input_dir $ctrl'_rc_all.txt' $input_dir $ct'.rep2'
	time Rscript ~/group/software/PKnorm/pknorm_scripts/fisher_pval.R $ct '.nbp_2r_bgadj.txt' $working_dir 323 $ct'.2r_nbp'
	tail -n+13274488 $ct'.2r_nbp.fisher_p.txt' > $ct'.2r_nbp.fisher_p.qtpcr.txt'
	head -13274487 $ct'.2r_nbp.fisher_p.txt' > $ct'.2r_nbp.fisher_p.200bp.txt'
done < data_info.txt

### normalize qtPCR & RNA-seq fisher's method: -log10(pvalue)
paste WT.2r_nbp.fisher_p.qtpcr.txt 0hr.2r_nbp.fisher_p.qtpcr.txt 4hr.2r_nbp.fisher_p.qtpcr.txt 6hr.2r_nbp.fisher_p.qtpcr.txt >  ALL.chipseq_sig.txt

### get qtPCR signal
tail -n+2 qt_pcr_sig.txt | sort -k3,3 -k4,4n | awk -F '\t' -v OFS='\t' '{print ($6+$7)/2}' > WT.qt_pcr_sig.txt
tail -n+2 qt_pcr_sig.txt | sort -k3,3 -k4,4n | awk -F '\t' -v OFS='\t' '{print ($8+$9)/2}' > 0hr.qt_pcr_sig.txt
tail -n+2 qt_pcr_sig.txt | sort -k3,3 -k4,4n | awk -F '\t' -v OFS='\t' '{print ($10+$11)/2}' > 4hr.qt_pcr_sig.txt
tail -n+2 qt_pcr_sig.txt | sort -k3,3 -k4,4n | awk -F '\t' -v OFS='\t' '{print ($12+$13)/2}' > 6hr.qt_pcr_sig.txt
paste WT.qt_pcr_sig.txt 0hr.qt_pcr_sig.txt 4hr.qt_pcr_sig.txt 6hr.qt_pcr_sig.txt >  ALL.qt_pcr_sig.txt

time Rscript get_qtPCR_norm.R
time Rscript get_label_bed.R

cp mm9_200.sort.label.bed IDEAS_run/
cd IDEAS_run/
qsub run_IDEAS.sh
