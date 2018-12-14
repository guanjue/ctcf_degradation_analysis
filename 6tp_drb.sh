###### get bed file of peak lists
### (1) get merged_peak list
cat *_CTCF.bed | sort -u | sort -k1,1 -k2,2n > all.bed

cat */*_CTCF*.bed | sort -u | sort -k1,1 -k2,2n > all.bed
#bedtools merge -i all.bed > allmerged0.bed
#cp all.bed allmerged.bed
cat all.bed | awk -F '\t' -v OFS='\t' '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1}' | sort -k1,1 -k2,2n > all.mid.bed
bedtools merge -i all.mid.bed -d 300 > all.merged.mid.bed
cat all.merged.mid.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-250, int(($2+$3)/2)+250}' > allmerged.bed

### get merge macs peak
cat allmerged.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' | sort -k4,4 > allmerged.sorted.bed
cat allmerged.sorted.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2)-500>=0) print $1, int(($2+$3)/2)-500, int(($2+$3)/2)+500, $1"_"$2"_"$3; else print $1, 0, int(($2+$3)/2)+500, $1"_"$2"_"$3}' | sort -k4,4 > allmerged.1kb.sorted.bed
cat allmerged.sorted.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2)-2500>=0) print $1, int(($2+$3)/2)-2500, int(($2+$3)/2)+2500, $1"_"$2"_"$3; else print $1, 0, int(($2+$3)/2)+2500, $1"_"$2"_"$3}' | sort -k4,4 > allmerged.5kb.sorted.bed
cat allmerged.sorted.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2)-5000>=0) print $1, int(($2+$3)/2)-5000, int(($2+$3)/2)+5000, $1"_"$2"_"$3; else print $1, 0, int(($2+$3)/2)+5000, $1"_"$2"_"$3}' | sort -k4,4 > allmerged.10kb.sorted.bed


### (3) get 500bp_win list
bedtools intersect -a /storage/home/gzx103/scratch/ctcf/6_time_points/mm9_500.bed -b allmerged.bed -v > mm9_500.nonpk.bed
Rscript /storage/home/gzx103/scratch/ctcf/6_time_points/get_random_500r.R
cat mm9_500.nonpk.50000.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' | sort -k4,4 > mm9_500.sorted.bed

cat mm9_500.sorted.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2)-500>=0) print $1, int(($2+$3)/2)-500, int(($2+$3)/2)+500, $1"_"$2"_"$3; else print $1, 0, int(($2+$3)/2)+500, $1"_"$2"_"$3}' | sort -k4,4 > mm9_500.1kb.sorted.bed
cat mm9_500.sorted.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2)-2500>=0) print $1, int(($2+$3)/2)-2500, int(($2+$3)/2)+2500, $1"_"$2"_"$3; else print $1, 0, int(($2+$3)/2)+2500, $1"_"$2"_"$3}' | sort -k4,4 > mm9_500.5kb.sorted.bed
cat mm9_500.sorted.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2)-5000>=0) print $1, int(($2+$3)/2)-5000, int(($2+$3)/2)+5000, $1"_"$2"_"$3; else print $1, 0, int(($2+$3)/2)+5000, $1"_"$2"_"$3}' | sort -k4,4 > mm9_500.10kb.sorted.bed


###### get reads count in sample and input for each peak lists
### (1) get the reads count for the merged_peak list
rm -r merge_peak_tab
mkdir merge_peak_tab
while read files
do
        id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
        input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
        echo $id
        echo $input_id
        ### get the foreground data
        time ~/group/software/ucsc/bigWigAverageOverBed $id'_CTCF.bw' allmerged.sorted.bed $id'.allmerged.sorted.tab'
        sort -k1,1 $id'.allmerged.sorted.tab' > $id'.allmerged.sorted.tab.tmp' && mv $id'.allmerged.sorted.tab.tmp' $id'.allmerged.sorted.tab'
        cat $id'.allmerged.sorted.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.allmerged.sorted.rc_per_bp.tab'
        mv $id'.allmerged.sorted.tab' merge_peak_tab
        ###### get background
        background_win=(1kb 5kb 10kb)
        #if [ ! -f $id'.allmerged.'$win'.sorted.rc_per_bp.tab' ]; then
                echo "First run"
                ### get window (1kb 5kb 10kb) input signal
                for win in "${background_win[@]}"
                do
                        echo $win
                        ### get background $win
                        time ~/group/software/ucsc/bigWigAverageOverBed $input_id'_Input.bw' 'allmerged.'$win'.sorted.bed' $id'.allmerged.'$win'.sorted.tab'
                        sort -k1,1 $id'.allmerged.'$win'.sorted.tab' > $id'.allmerged.'$win'.sorted.tab.tmp' && mv $id'.allmerged.'$win'.sorted.tab.tmp' $id'.allmerged.'$win'.sorted.tab'
                        cat $id'.allmerged.'$win'.sorted.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.allmerged.'$win'.sorted.rc_per_bp.tab'
                        mv $id'.allmerged.'$win'.sorted.tab' merge_peak_tab
                done
                ### get whole genome mean signal
                time ~/group/software/ucsc/bigWigAverageOverBed $input_id'_Input.bw' /storage/home/gzx103/scratch/ctcf/6_time_points/mm9.nochrM.wg.bed $input_id'.mm9.nochrM.wg.tab'
        #fi
        time Rscript /storage/home/gzx103/scratch/ctcf/6_time_points/get_local_bg_rc.R $input_id'.mm9.nochrM.wg.tab' $id'.allmerged.1kb.sorted.rc_per_bp.tab' $id'.allmerged.5kb.sorted.rc_per_bp.tab' $id'.allmerged.10kb.sorted.rc_per_bp.tab' $input_id'.allmerged.input_sig.macs.tab'
done < id_list_all.txt




### (3) get the reads count for the mm9_500 list
rm -r mm9_500_tab
mkdir mm9_500_tab
while read files
do
        id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $1}')
        input_id=$(echo "$files" | awk -F '\t' -v OFS='\t' '{print $2}')
        echo $id
        echo $input_id
        ### get the foreground data
        time ~/group/software/ucsc/bigWigAverageOverBed $id'_CTCF.bw' mm9_500.sorted.bed $id'.mm9_500.sorted.tab'
        sort -k1,1 $id'.mm9_500.sorted.tab' > $id'.mm9_500.sorted.tab.tmp' && mv $id'.mm9_500.sorted.tab.tmp' $id'.mm9_500.sorted.tab'
        cat $id'.mm9_500.sorted.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.mm9_500.sorted.rc_per_bp.tab'
        mv $id'.mm9_500.sorted.tab' mm9_500_tab
        ###### get background
        background_win=(1kb 5kb 10kb)
        #if [ ! -f $id'.mm9_500.'$win'.sorted.rc_per_bp.tab' ]; then
                echo "First run"
                ### get window (1kb 5kb 10kb) input signal
                for win in "${background_win[@]}"
                do
                        echo $win
                        ### get background $win
                        time ~/group/software/ucsc/bigWigAverageOverBed $input_id'_Input.bw' 'mm9_500.'$win'.sorted.bed' $id'.mm9_500.'$win'.sorted.tab'
                        sort -k1,1 $id'.mm9_500.'$win'.sorted.tab' > $id'.mm9_500.'$win'.sorted.tab.tmp' && mv $id'.mm9_500.'$win'.sorted.tab.tmp' $id'.mm9_500.'$win'.sorted.tab'
                        cat $id'.mm9_500.'$win'.sorted.tab' | awk -F '\t' -v OFS='\t' '{print $4/$2}' > $id'.mm9_500.'$win'.sorted.rc_per_bp.tab'
                        mv $id'.mm9_500.'$win'.sorted.tab' mm9_500_tab
                done
                ### get whole genome mean signal
                time ~/group/software/ucsc/bigWigAverageOverBed $input_id'_Input.bw' /storage/home/gzx103/scratch/ctcf/6_time_points/mm9.nochrM.wg.bed $input_id'.mm9.nochrM.wg.tab'
        #fi
        time Rscript /storage/home/gzx103/scratch/ctcf/6_time_points/get_local_bg_rc.R $input_id'.mm9.nochrM.wg.tab' $id'.mm9_500.1kb.sorted.rc_per_bp.tab' $id'.mm9_500.5kb.sorted.rc_per_bp.tab' $id'.mm9_500.10kb.sorted.rc_per_bp.tab' $input_id'.mm9_500.input_sig.macs.tab'
done < id_list_all.txt


