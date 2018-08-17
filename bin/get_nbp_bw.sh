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
	echo $ct
#	paste mm9_200.sort.bed $ct'.2r_nbp.fisher_p.200bp.normed.txt' > $ct'.2r_nbp.fisher_p.200bp.normed.bedgraph'
#	time /storage/home/gzx103/group/software/ucsc/bedGraphToBigWig $ct'.2r_nbp.fisher_p.200bp.normed.bedgraph' mm9.chrom_size.txt $ct'.S3norm.bw'
	paste mm9_200.sort.bed $ct'.2r_nbp.fisher_p.200bp.normed.txt' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,log($4+0.1)/log(2)}' > $ct'.2r_nbp.fisher_p.200bp.normed.log2.bedgraph'
	time /storage/home/gzx103/group/software/ucsc/bedGraphToBigWig $ct'.2r_nbp.fisher_p.200bp.normed.log2.bedgraph' mm9.chrom_size.txt $ct'.S3norm.log2.bw'
done < data_info.txt

