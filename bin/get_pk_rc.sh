#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc
module load python/2.7.14-anaconda5.0.1

cd /storage/home/gzx103/scratch/ctcf

#tail -n+2 ~/group/projects/vision/200_noblack.windows | awk -F ' ' -v OFS='\t' '{print $1,$2,$3}' > mm10_200_noblack.bed
head -21 IDEAS_run/data/mm9.genome > mm9.nochrM.genome
bedtools makewindows -g mm9.nochrM.genome -w 200 > mm9_200.bed
sort -k1,1 -k2,2n mm9_200.bed > mm9_200.sort.bed
#tail -n+2 qt_pcr_sig.txt | cut -f3,4,5 | sort -k1,1 -k2,2n > qtPCR_region.sort.bed
#sort -k3,3 -k4,4n qt_pcr_sig.txt > qt_pcr_sig.sort.txt

for id in $(cat id_list.txt)
do
	echo $id
#	time samtools index $id'_mouse_verysensitive1.filtered.bam'
#	time bamCoverage --binSize 200 --bam $id'_mouse_verysensitive1.filtered.bam' --outFileFormat bedgraph -o $id'_mouse_verysensitive1.filtered.bedgraph'
#	time /storage/home/gzx103/group/software/ucsc/bedGraphToBigWig $id'_mouse_verysensitive1.filtered.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $id'_mouse_verysensitive1.filtered.bw'

#	time bedtools map -a 201803106_AllCtCFPeaks_DiffBind.bed -b $id'_mouse_verysensitive1.filtered.bedgraph' -c 4 -o mean > $id'_peak_rc.txt'
	time bedtools map -a mm9_200.sort.bed -b $id'_mouse_verysensitive1.filtered.bedgraph' -c 4 -o mean -null 0 > $id'_200bp_rc.txt'
	time bedtools map -a qtPCR_region.sort.bed -b $id'_mouse_verysensitive1.filtered.bedgraph' -c 4 -o mean -null 0 > $id'_qtpcr_rc.txt'
done


### get qtPCR signal
tail -n+2 qt_pcr_sig.txt | sort -k3,3 -k4,4n | awk -F '\t' -v OFS='\t' '{print ($6+$7)/2}' > WT.qt_pcr_sig.txt
tail -n+2 qt_pcr_sig.txt | sort -k3,3 -k4,4n | awk -F '\t' -v OFS='\t' '{print ($8+$9)/2}' > 0hr.qt_pcr_sig.txt
tail -n+2 qt_pcr_sig.txt | sort -k3,3 -k4,4n | awk -F '\t' -v OFS='\t' '{print ($10+$11)/2}' > 4hr.qt_pcr_sig.txt
tail -n+2 qt_pcr_sig.txt | sort -k3,3 -k4,4n | awk -F '\t' -v OFS='\t' '{print ($12+$13)/2}' > 6hr.qt_pcr_sig.txt


qsub get_nbp.sh
