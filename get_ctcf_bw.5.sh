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
echo 'get read counts Fold-change with MACS BG'
###### get read counts Fold-change
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IP.bw' mm9_50.bin.bed $id'.mm9_50.sorted.tab'
	sort -k1,1 $id'.mm9_50.sorted.tab' > $id'.mm9_50.sorted.tab.tmp' && mv $id'.mm9_50.sorted.tab.tmp' $id'.mm9_50.sorted.tab'
	cat $id'.mm9_50.sorted.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.mm9_50.sorted.signal.1.tab'
	### get MACS input
	background_win=(1kb 5kb 10kb)
	echo "First run"
	### get window (1kb 5kb 10kb) input signal
	for win in "${background_win[@]}"
	do
		echo $win
		### get background $win
		time ~/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' 'mm9_50.bin.'$win'.bed' $id'.mm9_50.sorted.'$win'.tab'
		sort -k1,1 $id'.mm9_50.sorted.'$win'.tab' > $id'.mm9_50.sorted.'$win'.tab.tmp' && mv $id'.mm9_50.sorted.'$win'.tab.tmp' $id'.mm9_50.sorted.'$win'.tab'
		cat $id'.mm9_50.sorted.'$win'.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.mm9_50.sorted.signal.'$win'.tab'
	done
	### get whole genome mean signal
	time ~/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' mm9.nochrM.wg.bed $id'.mm9.nochrM.wg.tab'
	time Rscript get_local_bg_rc.R $id'.mm9.nochrM.wg.tab' $id'.mm9_50.sorted.signal.1kb.tab' $id'.mm9_50.sorted.signal.5kb.tab' $id'.mm9_50.sorted.signal.10kb.tab' $id'.mm9_50.sorted.signal.macsBG.tab'

	paste $id'.mm9_50.sorted.signal.1.tab' $id'.mm9_50.sorted.signal.macsBG.tab' | awk -F '\t' -v OFS='\t' '{if ($1/$2>=1) print $1/$2; else print 1}' > $id'.mm9_50.sorted.signal.fc.tab'
	rm $id'.mm9_50.sorted.signal.1.tab'
	rm $id'.mm9_50.sorted.signal.1kb.tab' $id'.mm9_50.sorted.signal.5kb.tab' $id'.mm9_50.sorted.signal.10kb.tab'
	rm $id'.mm9_50.sorted.1kb.tab' $id'.mm9_50.sorted.5kb.tab' $id'.mm9_50.sorted.10kb.tab'
	rm $id'.mm9.nochrM.wg.tab' $id'.mm9_50.sorted.signal.macsBG.tab'
	mv $id'.mm9_50.sorted.tab' bin_tab
	mv $id'.mm9_50.sorted.signal.fc.tab' bin_tab
done < bam_bw_bedgraph/file_list.all.5.txt




