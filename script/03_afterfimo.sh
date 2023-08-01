#!/bin/bash

# Processing fimo.gff file
#login fimo
awk '{print $1,$4,$5,$7,$9,$10}' $1|awk -F':' '{print $2,$3}'|awk -F';' '{print $1,$4,$3}'|sed 's/\([0-9]\+\)-\([0-9]\+\)/\1\t\2/g' |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$9}'|sed 's/pvalue=//g'  > fimo_processed.txt
#webdep fimo
#awk '{print $1,$4,$5,$7,$9,$10}' $1 |awk -F';' '{print $1,$4,$5}' |awk '{print $1,$2,$3,$4,$6,$7}' |awk -F':' '{print $2,$3}' | sed 's/pvalue=//g' |sed 's/qvalue=//g' |sed 's/\([0-9]\+\)-\([0-9]\+\)/\1\t\2/g'|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'  > fimo_processed.txt

head fimo_processed.txt 
# Header for the final output
echo -e "Strain\tMethylation\tRegion\tRRS\tRRE\tMethsite\tFraction\tdistance\tstart\tend\tstrand\tlocustag\tgene name\tdescription\tTF binding start\tTF binding end\tTF binding strand\tFIMO qvalue\tFIMO pvalue" > final_output.txt
awk 'NR==FNR{fimo[$1]=$0; next} ($4 in fimo) || ($5 in fimo){print $0 "\t" fimo[$4] "\t" fimo[$5]}' fimo_processed.txt $2_methylationGene.txt  > final_output_temp.txt
#awk 'NR==FNR{fimo[$1]=$0; next} ($4 in fimo) || ($5 in fimo){print $0 "\t" fimo[$4] "\t" fimo[$5]}' fimo_processed.txt $2_methylationGene.txt  > final_output_temp.txt

head final_output_temp.txt
awk  -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15+$17"\t"$15+$18"\t"$19"\t"$20"\t"$21}' final_output_temp.txt >>final_output.txt

awk -F'\t' '!($17 ~ /^[0-9]+$/)' final_output.txt >$2_${3}.txt
# Clean up
rm fimo_processed.txt final_output_temp.txt final_output.txt

