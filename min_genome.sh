paste mm9_200.sort.bed mm9.200bp.signal.mat.txt | awk -F '\t' -v OFS='\t' '{if (($4>2) || ($5>2) || ($6>2) || ($7>2)) print $0}'> mm9_200.sort.mini.txt

cut -f1,2,3 mm9_200.sort.mini.txt > mm9_200.sort.mini.bed

time get_label_bed.R

cut -f4 mm9_200.sort.mini.txt > WT.mm9_200.sort.mini.txt
cut -f5 mm9_200.sort.mini.txt > 0hr.mm9_200.sort.mini.txt
cut -f6 mm9_200.sort.mini.txt > 4hr.mm9_200.sort.mini.txt
cut -f7 mm9_200.sort.mini.txt > 6hr.mm9_200.sort.mini.txt

mv *.mm9_200.sort.mini.txt input_files/



