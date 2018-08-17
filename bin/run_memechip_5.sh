#PBS -l nodes=1:ppn=10
#PBS -l walltime=50:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc
module load python/2.7.14-anaconda5.0.1

cd /storage/home/gzx103/scratch/ctcf


while read LINE
do
	index=$(echo "$LINE" | awk -F '\t' '{print $1}')
	var1=$(echo "$LINE" | awk -F '\t' '{print $1}' | awk -F '_' '{print $1}')
	var2=$(echo "$LINE" | awk -F '\t' '{print $1}' | awk -F '_' -v OFS='_' '{print $2,$3,$4}')
	echo $index
	#echo $var1
	#echo $var2
	line_num=$(wc -l < $index'.bed')
	echo $line_num
	if (( line_num > 5000 )); then
		echo $index'.bed' 
		cat $index'.bed' | python -c "import random, sys; random.seed(2018); print ''.join(random.sample(sys.stdin.readlines(), 5000))," | awk -F '\t' -v OFS='\t' '{print $1,$2,$3, l2"_"l1";shared"}' > $index'.s.bed'
		bedtools getfasta -fi /storage/home/gzx103/group/genome/mm9/mm9.fasta -bed $index'.s.bed' > $index'.s.fa'
	else
		echo 'use all'
		cat $index'.bed' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3, l2"_"l1";shared"}' > $index'.s.bed'
		bedtools getfasta -fi /storage/home/gzx103/group/genome/mm9/mm9.fasta -bed $index'.s.bed' > $index'.s.fa'
	fi 
done < index_count_enrich_ic_noXXXX.txt

