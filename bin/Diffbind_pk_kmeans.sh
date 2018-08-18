time Rscript cluster_index_set_pk_kmeans.R

paste 201803106_AllCtCFPeaks_DiffBind.nochrUn_random.bed kmeans_id_index_vec.txt 201803106_AllCtCFPeaks_DiffBind.S3norm.nochrUn_random.txt > 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.txt

cat 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4, $5,$6,$7,$8, log($5)+log($6)+log($7)+log($8)}' | sort -k4,4r -k9,9rn > 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.sort.txt

cut -f1,2,3,4,5,6,7,8 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.sort.txt > 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.sort.meansigsort.txt

time Rscript /storage/home/gzx103/group/software/ctcf_degradation_analysis/bin/box_trend.R

#cat 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.meansigsort.txt | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3, $4, $5, $6, $7, $8}' > 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.meansigsort.foric.txt

#time Rscript ~/group/software/snapshot/bin/plot_rect_sig.R 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.meansigsort.txt 201803106_AllCtCFPeaks_DiffBind.S3norm.index.sort.png tp_list.txt 5 red white gray T 0.1

cut -f1,2,3 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.sort.txt > 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.sort.meansigsort.bed

while read LINE
do
	index=$(echo "$LINE" | awk -F '\t' '{print $1}')
	echo $index
	cat 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.sort.txt | awk -F '\t' -v OFS='\t' -v index_tar="$index" '{if ($4==index_tar) print $0}' > '201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.'$index'.txt'
	cut -f1,2,3 '201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.'$index'.txt' > '201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.'$index'.bed'
	time computeMatrix scale-regions -S WT.S3norm.log2.bw 0hr.S3norm.log2.bw 4hr.S3norm.log2.bw 6hr.S3norm.log2.bw -R '201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.'$index'.bed' --beforeRegionStartLength 10000 --regionBodyLength 1000 --afterRegionStartLength 10000 -o '201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.'$index'.gz' --binSize 200 --numberOfProcessors 5 --sortRegions keep --missingDataAsZero --averageTypeBins mean
	time plotHeatmap --colorList 'white, red' -m '201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.'$index'.gz' -out '201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.'$index'.png' --sortRegions no --zMax 9 --zMin -2 --yMin -2 --yMax 9
done < kmeans_id_index_count.txt


ls -r *201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.*.bed > bed_list_kmeans.txt

for file in $(cat bed_list_kmeans.txt)
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
done < kmeans_id_index_count.txt


time computeMatrix scale-regions -S WT.S3norm.log2.bw 0hr.S3norm.log2.bw 4hr.S3norm.log2.bw 6hr.S3norm.log2.bw -R $foo --beforeRegionStartLength 10000 --regionBodyLength 1000 --afterRegionStartLength 10000 -o 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.enriched_index.gz --binSize 200 --numberOfProcessors 5 --sortRegions keep --missingDataAsZero --averageTypeBins mean

time plotHeatmap --colorList 'white, red' -m 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.enriched_index.gz -out 201803106_AllCtCFPeaks_DiffBind.S3norm.kmeans.enriched_index.pdf --sortRegions no --zMax 9 --zMin -2 --yMin -2 --yMax 9 --startLabel s --endLabel e --heatmapHeight 60







