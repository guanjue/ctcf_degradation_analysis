############################################################
############################################################
############################################################
rm -r qtPCR_regions_tab
mkdir qtPCR_regions_tab
rm -r qtPCR_regions_signal_tab
mkdir qtPCR_regions_signal_tab

####################################
echo 'get read counts Fold-change'
###### get read counts Fold-change
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IP.bw' qtPCR_regions.sorted.5T7D.bed $id'.qtPCR_regions.sorted.5T7D.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.5T7D.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.5T7D.tab'
	cat $id'.qtPCR_regions.sorted.5T7D.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.signal.5T7D1.tab'
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' qtPCR_regions.sorted.5T7D.bed $id'.qtPCR_regions.sorted.5T7D.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.5T7D.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.5T7D.tab'
	cat $id'.qtPCR_regions.sorted.5T7D.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.signal.5T7D2.tab'
	paste $id'.qtPCR_regions.sorted.signal.5T7D1.tab' $id'.qtPCR_regions.sorted.signal.5T7D2.tab' | awk -F '\t' -v OFS='\t' '{print ($1+1)/($2+1)}' > $id'.qtPCR_regions.sorted.signal.5T7D.fc0.tab'
	rm $id'.qtPCR_regions.sorted.signal.5T7D1.tab' $id'.qtPCR_regions.sorted.signal.5T7D2.tab'
	mv $id'.qtPCR_regions.sorted.5T7D.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.5T7D.fc0.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt

### (2) get the reads count for the qtPCR_regions list
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IP.bw' qtPCR_regions.sorted.6TP.bed $id'.qtPCR_regions.sorted.6TP.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.6TP.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.6TP.tab'
	cat $id'.qtPCR_regions.sorted.6TP.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.signal.6TP1.tab'
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' qtPCR_regions.sorted.6TP.bed $id'.qtPCR_regions.sorted.6TP.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.6TP.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.6TP.tab'
	cat $id'.qtPCR_regions.sorted.6TP.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.qtPCR_regions.sorted.signal.6TP2.tab'
	paste $id'.qtPCR_regions.sorted.signal.6TP1.tab' $id'.qtPCR_regions.sorted.signal.6TP2.tab' | awk -F '\t' -v OFS='\t' '{print ($1+1)/($2+1)}' > $id'.qtPCR_regions.sorted.signal.6TP.fc0.tab'
	rm $id'.qtPCR_regions.sorted.signal.6TP1.tab' $id'.qtPCR_regions.sorted.signal.6TP2.tab'
	mv $id'.qtPCR_regions.sorted.6TP.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.6TP.fc0.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt


####################################
echo 'get read counts log2 Fold-change'
###### get read counts log2 Fold-change
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IPlogfcCTRL_SES.pos.bw' qtPCR_regions.sorted.5T7D.bed $id'.qtPCR_regions.sorted.5T7D.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.5T7D.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.5T7D.tab'
	cat $id'.qtPCR_regions.sorted.5T7D.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.signal.5T7D.log2fc.tab'
	mv $id'.qtPCR_regions.sorted.5T7D.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.5T7D.log2fc.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt

### (2) get the reads count for the qtPCR_regions list
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IPlogfcCTRL_SES.pos.bw' qtPCR_regions.sorted.6TP.bed $id'.qtPCR_regions.sorted.6TP.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.6TP.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.6TP.tab'
	cat $id'.qtPCR_regions.sorted.6TP.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.signal.6TP.log2fc.tab'
	mv $id'.qtPCR_regions.sorted.6TP.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.6TP.log2fc.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt


####################################
echo 'get read counts subtract'
###### get read counts subtract
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IPsubCTRL_SES.pos.bw' qtPCR_regions.sorted.5T7D.bed $id'.qtPCR_regions.sorted.5T7D.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.5T7D.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.5T7D.tab'
	cat $id'.qtPCR_regions.sorted.5T7D.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.signal.5T7D.sub.tab'
	mv $id'.qtPCR_regions.sorted.5T7D.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.5T7D.sub.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt

### (2) get the reads count for the qtPCR_regions list
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IPsubCTRL_SES.pos.bw' qtPCR_regions.sorted.6TP.bed $id'.qtPCR_regions.sorted.6TP.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.6TP.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.6TP.tab'
	cat $id'.qtPCR_regions.sorted.6TP.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.signal.6TP.sub.tab'
	mv $id'.qtPCR_regions.sorted.6TP.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.6TP.sub.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt


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
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IP.bw' qtPCR_regions.sorted.5T7D.bed $id'.qtPCR_regions.sorted.5T7D.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.5T7D.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.5T7D.tab'
	cat $id'.qtPCR_regions.sorted.5T7D.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.qtPCR_regions.sorted.signal.5T7D1.tab'
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' qtPCR_regions.sorted.5T7D.bed $id'.qtPCR_regions.sorted.5T7D.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.5T7D.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.5T7D.tab'
	cat $id'.qtPCR_regions.sorted.5T7D.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.qtPCR_regions.sorted.signal.5T7D2.tab'

	background_win=(1kb 5kb 10kb)
	echo "First run"
	### get window (1kb 5kb 10kb) input signal
	for win in "${background_win[@]}"
	do
		echo $win
		### get background $win
		time ~/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' 'qtPCR_regions.sorted.5T7D.'$win'.bed' $id'.qtPCR_regions.sorted.5T7D.'$win'.tab'
		sort -k1,1 $id'.qtPCR_regions.sorted.5T7D.'$win'.tab' > $id'.qtPCR_regions.sorted.5T7D.'$win'.tab.tmp' && mv $id'.qtPCR_regions.sorted.5T7D.'$win'.tab.tmp' $id'.qtPCR_regions.sorted.5T7D.'$win'.tab'
		cat $id'.qtPCR_regions.sorted.5T7D.'$win'.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.qtPCR_regions.sorted.signal.5T7D.'$win'.tab'
	done
	### get whole genome mean signal
	time ~/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' mm9.nochrM.wg.bed $id'.mm9.nochrM.wg.tab'
	time Rscript get_local_bg_rc.R $id'.mm9.nochrM.wg.tab' $id'.qtPCR_regions.sorted.signal.5T7D.1kb.tab' $id'.qtPCR_regions.sorted.signal.5T7D.5kb.tab' $id'.qtPCR_regions.sorted.signal.5T7D.10kb.tab' $id'.qtPCR_regions.sorted.signal.5T7D.macsBG.tab'

	paste $id'.qtPCR_regions.sorted.signal.5T7D1.tab' $id'.qtPCR_regions.sorted.signal.5T7D.macsBG.tab' | awk -F '\t' -v OFS='\t' '{if ($1/$2>=1) print $1/$2; else print 1}' > $id'.qtPCR_regions.sorted.signal.5T7D.fc.tab'
	rm $id'.qtPCR_regions.sorted.signal.5T7D1.tab' $id'.qtPCR_regions.sorted.signal.5T7D2.tab' 
	rm $id'.qtPCR_regions.sorted.signal.5T7D.1kb.tab' $id'.qtPCR_regions.sorted.signal.5T7D.5kb.tab' $id'.qtPCR_regions.sorted.signal.5T7D.10kb.tab'
	rm $id'.qtPCR_regions.sorted.5T7D.1kb.tab' $id'.qtPCR_regions.sorted.5T7D.5kb.tab' $id'.qtPCR_regions.sorted.5T7D.10kb.tab'
	rm $id'.mm9.nochrM.wg.tab' $id'.qtPCR_regions.sorted.signal.5T7D.macsBG.tab'
	mv $id'.qtPCR_regions.sorted.5T7D.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.5T7D.fc.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt


### (2) get the reads count for the qtPCR_regions list
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IP.bw' qtPCR_regions.sorted.6TP.bed $id'.qtPCR_regions.sorted.6TP.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.6TP.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.6TP.tab'
	cat $id'.qtPCR_regions.sorted.6TP.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.qtPCR_regions.sorted.signal.6TP1.tab'
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' qtPCR_regions.sorted.6TP.bed $id'.qtPCR_regions.sorted.6TP.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.6TP.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.6TP.tab'
	cat $id'.qtPCR_regions.sorted.6TP.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.qtPCR_regions.sorted.signal.6TP2.tab'

	background_win=(1kb 5kb 10kb)
	echo "First run"
	### get window (1kb 5kb 10kb) input signal
	for win in "${background_win[@]}"
	do
		echo $win
		### get background $win
		time ~/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' 'qtPCR_regions.sorted.6TP.'$win'.bed' $id'.qtPCR_regions.sorted.6TP.'$win'.tab'
		sort -k1,1 $id'.qtPCR_regions.sorted.6TP.'$win'.tab' > $id'.qtPCR_regions.sorted.6TP.'$win'.tab.tmp' && mv $id'.qtPCR_regions.sorted.6TP.'$win'.tab.tmp' $id'.qtPCR_regions.sorted.6TP.'$win'.tab'
		cat $id'.qtPCR_regions.sorted.6TP.'$win'.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.qtPCR_regions.sorted.signal.6TP.'$win'.tab'
	done
	### get whole genome mean signal
	time ~/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' mm9.nochrM.wg.bed $id'.mm9.nochrM.wg.tab'
	time Rscript get_local_bg_rc.R $id'.mm9.nochrM.wg.tab' $id'.qtPCR_regions.sorted.signal.6TP.1kb.tab' $id'.qtPCR_regions.sorted.signal.6TP.5kb.tab' $id'.qtPCR_regions.sorted.signal.6TP.10kb.tab' $id'.qtPCR_regions.sorted.signal.6TP.macsBG.tab'

	paste $id'.qtPCR_regions.sorted.signal.6TP1.tab' $id'.qtPCR_regions.sorted.signal.6TP.macsBG.tab' | awk -F '\t' -v OFS='\t' '{if ($1/$2>=1) print $1/$2; else print 1}' > $id'.qtPCR_regions.sorted.signal.6TP.fc.tab'
	rm $id'.qtPCR_regions.sorted.signal.6TP1.tab' $id'.qtPCR_regions.sorted.signal.6TP2.tab' 
	rm $id'.qtPCR_regions.sorted.signal.6TP.1kb.tab' $id'.qtPCR_regions.sorted.signal.6TP.5kb.tab' $id'.qtPCR_regions.sorted.signal.6TP.10kb.tab'
	rm $id'.qtPCR_regions.sorted.6TP.1kb.tab' $id'.qtPCR_regions.sorted.6TP.5kb.tab' $id'.qtPCR_regions.sorted.6TP.10kb.tab'
	rm $id'.mm9.nochrM.wg.tab' $id'.qtPCR_regions.sorted.signal.6TP.macsBG.tab'
	mv $id'.qtPCR_regions.sorted.6TP.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.6TP.fc.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt


####################################
echo 'get read counts RAW'
###### get read counts RAW
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IP.bw' qtPCR_regions.sorted.5T7D.bed $id'.qtPCR_regions.sorted.5T7D.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.5T7D.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.5T7D.tab'
	cat $id'.qtPCR_regions.sorted.5T7D.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.signal.5T7D.raw.tab'
	mv $id'.qtPCR_regions.sorted.5T7D.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.5T7D.raw.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt

### (2) get the reads count for the qtPCR_regions list
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IP.bw' qtPCR_regions.sorted.6TP.bed $id'.qtPCR_regions.sorted.6TP.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.6TP.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.6TP.tab'
	cat $id'.qtPCR_regions.sorted.6TP.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.signal.6TP.raw.tab'
	mv $id'.qtPCR_regions.sorted.6TP.tab' qtPCR_regions_tab
	mv $id'.qtPCR_regions.sorted.signal.6TP.raw.tab' qtPCR_regions_signal_tab
done < bam_bw_bedgraph/file_list.all.txt




