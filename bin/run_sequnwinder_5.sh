#PBS -l nodes=1:ppn=10
#PBS -l walltime=50:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc
module load python/2.7.14-anaconda5.0.1

cd /storage/home/gzx103/scratch/ctcf


rm all_cluster_s_5_ic.txt

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
		cat $index'.bed' | python -c "import random, sys; random.seed(2018); print ''.join(random.sample(sys.stdin.readlines(), 5000))," | awk -F '\t' -v OFS='\t' -v l1="$var1" -v l2="$var2" '{print $1":"int(($2+$3)/2), l2"_"l1";shared"}' >> all_cluster_s_5_ic.txt
	else
		echo 'use all'
		cat $index'.bed' | awk -F '\t' -v OFS='\t' -v l1="$var1" -v l2="$var2" '{print $1":"int(($2+$3)/2), l2"_"l1";shared"}' >> all_cluster_s_5_ic.txt
	fi 
done < index_count_enrich_ic_noXXXX.txt


time java -Xmx20G -jar /storage/home/gzx103/group/software/sequnwinder/sequnwinder_v0.1.3.jar --out ctcf_auxin_5_ic --threads 10 --debug --memepath /storage/home/gzx103/group/software/meme/bin --geninfo mm9.info --seq /storage/home/gzx103/group/genome/mm9/ --genregs all_cluster_s_5.txt --win 500 --mink 4 --maxk 5 --r 10 --x 3 --maxscanlen 15



### min version

rm all_cluster_s_5_ic_mini.txt

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
if (( line_num > 500 )); then
echo $index'.bed' 
cat $index'.bed' | python -c "import random, sys; random.seed(2018); print ''.join(random.sample(sys.stdin.readlines(), 500))," | awk -F '\t' -v OFS='\t' -v l1="$var1" -v l2="$var2" '{print $1":"int(($2+$3)/2), l2"_"l1";shared"}' >> all_cluster_s_5_ic_mini.txt
else
echo 'use all'
cat $index'.bed' | awk -F '\t' -v OFS='\t' -v l1="$var1" -v l2="$var2" '{print $1":"int(($2+$3)/2), l2"_"l1";shared"}' >> all_cluster_s_5_ic_mini.txt
fi 
done < index_count_enrich_ic_noXXXX.txt


time java -Xmx20G -jar /storage/home/gzx103/group/software/sequnwinder/sequnwinder_v0.1.3.jar --out ctcf_auxin_5_ic_mini --threads 10 --debug --memepath /storage/home/gzx103/group/software/meme/bin --geninfo mm9.info --seq /storage/home/gzx103/group/genome/mm9/ --genregs all_cluster_s_5_ic_mini.txt --win 500 --mink 4 --maxk 5 --r 10 --x 3 --maxscanlen 15


