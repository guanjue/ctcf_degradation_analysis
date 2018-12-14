#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

source ~/.bash_profile

cd /storage/home/gzx103/group/projects/ctcf_auxin/

####################################
echo 'Get normalized whole genome bins'
###### Get normalized whole genome bins
time Rscript LM_s3_6TP_wgbin.R


####################################
echo 'get read counts Fold-change with MACS BG'
###### get read counts Fold-change
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### sort peak
	sort -k1,1 -k2,2n 'bin_tab/'$id'.mm9_50.sorted.signal.fc.6TP.bedgraph' > 'bin_tab/'$id'.mm9_50.sorted.signal.fc.sorted.6TP.bedgraph'
	### get bw files
	bedGraphToBigWig 'bin_tab/'$id'.mm9_50.sorted.signal.fc.sorted.6TP.bedgraph' ~/group/genome/mm9/mm9.1to19_X.genome 'bin_tab/'$id'.mm9_50.6TP.fcMACS.bw'
done < bam_bw_bedgraph/file_list.all.txt


