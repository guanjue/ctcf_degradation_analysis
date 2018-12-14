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


############################################################
#######(((1)))###### get bed file of peak lists
### (1) get merged_peak list
cat peaks/*_CTCF*.bed | sort -u | sort -k1,1 -k2,2n > all.bed
### get mid point
cat all.bed | awk -F '\t' -v OFS='\t' '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1}' | sort -k1,1 -k2,2n > all.mid.bed
### merge mid point if distance <= 250-bp
bedtools merge -i all.mid.bed -d 250 > all.merged.mid.bed
### expand mid point to peak (+- 250-bp)
cat all.merged.mid.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-250, int(($2+$3)/2)+250, $1"_"int(($2+$3)/2)-250"_"int(($2+$3)/2)+250}' > allmerged.bed
sort -k1,1 -k2,2n allmerged.bed > allmerged.sorted.bed


###### get reads count in sample and input for each peak lists
### (2) get the reads count for the merged_peak list
rm -r merge_peak_tab
mkdir merge_peak_tab
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	#input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	#echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IPsubCTRL_SES.pos.bw' allmerged.sorted.bed $id'.allmerged.sorted.tab'
	sort -k1,1 $id'.allmerged.sorted.tab' > $id'.allmerged.sorted.tab.tmp' && mv $id'.allmerged.sorted.tab.tmp' $id'.allmerged.sorted.tab'
	cat $id'.allmerged.sorted.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.allmerged.sorted.rc_per_bp.tab'
	mv $id'.allmerged.sorted.tab' merge_peak_tab
done < bam_bw_bedgraph/file_list.all.txt


############################################################
#######(((2)))###### get bed file of qPCR lists
### (1) extract the bed file of the qPCR table
tail -n+2 20181206_DRB_Triptolide_ChIPqPCR.trimmed.csv | awk -F ',' -v OFS='\t' '{print $1, $2, $3, $4}' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3, $4}' | sort -k4,4 > qtPCR_regions.sorted.bed.txt
cut -f1,2,3,4 qtPCR_regions.sorted.bed.txt > qtPCR_regions.sorted.bed

### (2) get the reads count for the qtPCR_regions list
rm -r qtPCR_regions_tab
mkdir qtPCR_regions_tab
while read files
do
	id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
	input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
	echo $id
	echo $input_id
	### get the foreground data
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IPsubCTRL_SES.pos.bw' qtPCR_regions.sorted.bed $id'.qtPCR_regions.sorted.tab'
	sort -k1,1 $id'.qtPCR_regions.sorted.tab' > $id'.qtPCR_regions.sorted.tab.tmp' && mv $id'.qtPCR_regions.sorted.tab.tmp' $id'.qtPCR_regions.sorted.tab'
	cat $id'.qtPCR_regions.sorted.tab' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'.qtPCR_regions.sorted.rc_per_bp.tab'
	mv $id'.qtPCR_regions.sorted.tab' qtPCR_regions_tab
done < bam_bw_bedgraph/file_list.all.txt


###### mv mkdir rc_per_bp.tab
rm -r rc_per_bp_tab
mkdir rc_per_bp_tab
mv *.rc_per_bp.tab rc_per_bp_tab

 




