paste 201803106_AllCtCFPeaks_DiffBind.nochrUn_random.bed signal_mat_index_vec.txt 201803106_AllCtCFPeaks_DiffBind.S3norm.nochrUn_random.txt > 201803106_AllCtCFPeaks_DiffBind.S3norm.index.txt

cat 201803106_AllCtCFPeaks_DiffBind.S3norm.index.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4, $5,$6,$7,$8, log($5)+log($6)+log($7)+log($8)}' | sort -k4,4r -k9,9rn > 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.txt

cut -f1,2,3,4,5,6,7,8 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.txt > 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.meansigsort.txt

cat 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.meansigsort.txt | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3, $4, $5, $6, $7, $8}' > 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.meansigsort.foric.txt

time Rscript ~/group/software/snapshot/bin/plot_rect_sig.R 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.meansigsort.txt 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.png tp_list.txt 5 red white gray T 0.1

cut -f1,2,3 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.txt > 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.meansigsort.bed

while read LINE
do
	index=$(echo "$LINE" | awk -F '\t' '{print $1}')
	echo $index
	cat 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.txt | awk -F '\t' -v OFS='\t' -v index_tar="$index" '{if ($4==index_tar) print $0}' > '201803106_AllCtCFPeaks_DiffBind.S3norm.index.'$index'.txt'
	cut -f1,2,3 '201803106_AllCtCFPeaks_DiffBind.S3norm.index.'$index'.txt' > '201803106_AllCtCFPeaks_DiffBind.S3norm.index.'$index'.bed'
	time computeMatrix scale-regions -S WT.S3norm.log2.bw 0hr.S3norm.log2.bw 4hr.S3norm.log2.bw 6hr.S3norm.log2.bw -R '201803106_AllCtCFPeaks_DiffBind.S3norm.index.'$index'.bed' --beforeRegionStartLength 10000 --regionBodyLength 1000 --afterRegionStartLength 10000 -o '201803106_AllCtCFPeaks_DiffBind.S3norm.index.'$index'.gz' --binSize 200 --numberOfProcessors max --sortRegions keep --missingDataAsZero --averageTypeBins mean
	time plotHeatmap --colorList 'white, red' -m '201803106_AllCtCFPeaks_DiffBind.S3norm.index.'$index'.gz' -out '201803106_AllCtCFPeaks_DiffBind.S3norm.index.'$index'.png' --sortRegions no --zMax 9 --zMin -2 --yMin -2 --yMax 9
done < index_count_enrich.txt


ls -r *201803106_AllCtCFPeaks_DiffBind.S3norm.index.*.bed > bed_list.txt

for file in $(cat bed_list.txt)
do
	index=$(echo "$file" | awk -F '.' '{print $4}')
	echo $index
	mv $file $index'.bed'
done

foo=""
while read LINE
do
	index=$(echo "$LINE" | awk -F '\t' '{print $1".bed"}')
	echo $index
	foo="$foo $index"
done < index_count_enrich.txt


time computeMatrix scale-regions -S WT.S3norm.log2.bw 0hr.S3norm.log2.bw 4hr.S3norm.log2.bw 6hr.S3norm.log2.bw -R $foo --beforeRegionStartLength 10000 --regionBodyLength 1000 --afterRegionStartLength 10000 -o 201803106_AllCtCFPeaks_DiffBind.S3norm.index.enriched_index.gz --binSize 200 --numberOfProcessors max --sortRegions keep --missingDataAsZero --averageTypeBins mean








