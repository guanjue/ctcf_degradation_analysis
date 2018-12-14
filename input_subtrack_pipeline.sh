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
cat all.merged.mid.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-250, int(($2+$3)/2)+250}' > allmerged.bed
sort -k1,1 -k2,2n allmerged.bed > allmerged.sorted.bed
cat allmerged.sorted.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-500, int(($2+$3)/2)+500, $1"_"$2"_"$3}' > allmerged.sorted.1kb.bed
cat allmerged.sorted.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-2500, int(($2+$3)/2)+2500, $1"_"$2"_"$3}' > allmerged.sorted.5kb.bed
cat allmerged.sorted.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-5000, int(($2+$3)/2)+5000, $1"_"$2"_"$3}' > allmerged.sorted.10kb.bed


###### get reads count in sample and input for each peak lists
### (2) get the reads count for the merged_peak list
rm -r merge_peak_tab
mkdir merge_peak_tab

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
	time /storage/home/gzx103/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'IP.bw' allmerged.sorted.bed $id'.allmerged.sorted.tab'
	sort -k1,1 $id'.allmerged.sorted.tab' > $id'.allmerged.sorted.tab.tmp' && mv $id'.allmerged.sorted.tab.tmp' $id'.allmerged.sorted.tab'
	cat $id'.allmerged.sorted.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.allmerged.sorted.signal.1.tab'
	### get MACS input
	background_win=(1kb 5kb 10kb)
	echo "First run"
	### get window (1kb 5kb 10kb) input signal
	for win in "${background_win[@]}"
	do
		echo $win
		### get background $win
		time ~/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' 'allmerged.sorted.'$win'.bed' $id'.allmerged.sorted.'$win'.tab'
		sort -k1,1 $id'.allmerged.sorted.'$win'.tab' > $id'.allmerged.sorted.'$win'.tab.tmp' && mv $id'.allmerged.sorted.'$win'.tab.tmp' $id'.allmerged.sorted.'$win'.tab'
		cat $id'.allmerged.sorted.'$win'.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.allmerged.sorted.signal.'$win'.tab'
	done
	### get whole genome mean signal
	time ~/group/software/ucsc/bigWigAverageOverBed 'bam_bw_bedgraph/'$id'CTRL.bw' mm9.nochrM.wg.bed $id'.mm9.nochrM.wg.tab'
	time Rscript get_local_bg_rc.R $id'.mm9.nochrM.wg.tab' $id'.allmerged.sorted.signal.1kb.tab' $id'.allmerged.sorted.signal.5kb.tab' $id'.allmerged.sorted.signal.10kb.tab' $id'.allmerged.sorted.signal.macsBG.tab'

	paste $id'.allmerged.sorted.signal.1.tab' $id'.allmerged.sorted.signal.macsBG.tab' | awk -F '\t' -v OFS='\t' '{print ($1)/($2)}' > $id'.allmerged.sorted.signal.fc.tab'
	rm $id'.allmerged.sorted.signal.1.tab'
	rm $id'.allmerged.sorted.signal.1kb.tab' $id'.allmerged.sorted.signal.5kb.tab' $id'.allmerged.sorted.signal.10kb.tab'
	rm $id'.allmerged.sorted.1kb.tab' $id'.allmerged.sorted.5kb.tab' $id'.allmerged.sorted.10kb.tab'
	rm $id'.mm9.nochrM.wg.tab' $id'.allmerged.sorted.signal.macsBG.tab'
	mv $id'.allmerged.sorted.tab' merge_peak_tab
	mv $id'.allmerged.sorted.signal.fc.tab' merge_peak_tab
done < bam_bw_bedgraph/file_list.all.txt





############################################################ ###### FOR 5T/7D ######
#######(((2)))###### get bed file of qPCR lists ###### FOR 5T/7D ######
### (1) extract the bed file of the qPCR table
tail -n+2 20181206_DRB_Triptolide_ChIPqPCR.trimmed.txt | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3, $4}' | sort -k4,4 > qtPCR_regions.sorted.5T7D.bed.txt
cut -f1,2,3,4 qtPCR_regions.sorted.5T7D.bed.txt > qtPCR_regions.sorted.5T7D.bed
cat qtPCR_regions.sorted.5T7D.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-500, int(($2+$3)/2)+500, $1"_"$2"_"$3}' > qtPCR_regions.sorted.5T7D.1kb.bed
cat qtPCR_regions.sorted.5T7D.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-2500, int(($2+$3)/2)+2500, $1"_"$2"_"$3}' > qtPCR_regions.sorted.5T7D.5kb.bed
cat qtPCR_regions.sorted.5T7D.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-5000, int(($2+$3)/2)+5000, $1"_"$2"_"$3}' > qtPCR_regions.sorted.5T7D.10kb.bed


############################################################
###### get qtPCR signal array
head -1 20181206_DRB_Triptolide_ChIPqPCR.trimmed.txt | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3, $4, $0}' > 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt
tail -n+2 20181206_DRB_Triptolide_ChIPqPCR.trimmed.txt | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3, $4, $0}' | sort -k4,4 >> 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt

tail -n+2 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($11+$12)/2*1000/($3-$2)}' > qtPCR_sig.5T7D.0A.txt
tail -n+2 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($13+$14+$15)/3*1000/($3-$2)}' > qtPCR_sig.5T7D.5T.txt
tail -n+2 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($16+$17)/2*1000/($3-$2)}' > qtPCR_sig.5T7D.7D.txt
tail -n+2 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($18+$19)/2*1000/($3-$2)}' > qtPCR_sig.5T7D.4A.txt
tail -n+2 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($20+$21+$22+$23)/4*1000/($3-$2)}' > qtPCR_sig.5T7D.5T4A.txt
tail -n+2 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($24+$25+$26+$27)/4*1000/($3-$2)}' > qtPCR_sig.5T7D.7D4A.txt
tail -n+2 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print $10}' > qtPCR_sig.5T7D.ids.txt
tail -n+2 20181206_DRB_Triptolide.trimmed.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print $3-$2}' > qtPCR_sig.5T7D.len.txt


############################################################ ###### FOR 6TP ######
#######(((3)))###### get bed file of qPCR lists ###### FOR 6TP ######
### (1) extract the bed file of the qPCR table
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.txt | awk -F '\t' -v OFS='\t' '{print $2, $3, $4, $2"_"$3"_"$4, $5}' | sort -k4,4 > qtPCR_regions.sorted.6TP.bed.txt
cut -f1,2,3,4 qtPCR_regions.sorted.6TP.bed.txt > qtPCR_regions.sorted.6TP.bed
cat qtPCR_regions.sorted.6TP.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-500, int(($2+$3)/2)+500, $1"_"$2"_"$3}' > qtPCR_regions.sorted.6TP.1kb.bed
cat qtPCR_regions.sorted.6TP.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-2500, int(($2+$3)/2)+2500, $1"_"$2"_"$3}' > qtPCR_regions.sorted.6TP.5kb.bed
cat qtPCR_regions.sorted.6TP.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-5000, int(($2+$3)/2)+5000, $1"_"$2"_"$3}' > qtPCR_regions.sorted.6TP.10kb.bed

############################################################
###### get qtPCR signal array
head -1 20180910_CTCF_ChIPqPCR_7pts.txt | awk -F '\t' -v OFS='\t' '{print $2, $3, $4, $2"_"$3"_"$4, $5, $0}' > 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.txt | awk -F '\t' -v OFS='\t' '{print $2, $3, $4, $2"_"$3"_"$4, $5, $0}' | sort -k4,4 >> 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt

tail -n+2 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($13+$14)/2*1000/($3-$2)}' > qtPCR_sig.6TP.0A.txt
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($15+$16)/2*1000/($3-$2)}' > qtPCR_sig.6TP.4A.txt
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($17+$18)/2*1000/($3-$2)}' > qtPCR_sig.6TP.6A.txt
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($19+$20)/2*1000/($3-$2)}' > qtPCR_sig.6TP.12A.txt
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($21+$22)/2*1000/($3-$2)}' > qtPCR_sig.6TP.18A.txt
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print ($23+$24)/2*1000/($3-$2)}' > qtPCR_sig.6TP.24A.txt
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print $6}' > qtPCR_sig.6TP.ids.txt
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.qtPCR_regions.sorted.bed.allinfo.txt | awk -F '\t' -v OFS='\t' '{print $3-$2}' > qtPCR_sig.6TP.len.txt


############################################################
###### get whole genome bed
cat mm9.nochrM.genome | awk -F '\t' -v OFS='\t' '{print $1, 0, $2, $1"_0_"$2}' | sort -k4,4 > mm9.nochrM.wg.bed


############################################################
###### get qPCR region ChIP-seq signal
time bash get_qPCR_signal.sh


############################################################
###### evaluate qPCR VS ChIP-seq correlation
time Rscript get_qPCR_norm.6TP.R
time Rscript get_qPCR_norm.5T7D.R


############################################################
###### get the final result peak (peak name sorted)
cut -f1 merge_peak_tab/1579_0A_CTCF_mm9_duprm.bam.allmerged.sorted.tab | awk -F '_' -v OFS='\t' '{print $1, $2, $3}' > allmerged.peakname_sorted.bed
cut -f1 merge_peak_tab/1579_0A_CTCF_mm9_duprm.bam.allmerged.sorted.tab > allmerged.peakname_sorted.bed.txt

############################################################
###### normalize merge_peak_tab signal












