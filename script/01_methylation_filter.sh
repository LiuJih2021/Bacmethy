#!/usr/bin/env bash
#liujihong
coverage=$6
#echo $*
#echo $#
#echo $coverage
cat $1 | tr "," "\t" > test.tsv 
grep "$3" test.tsv > "$3".tsv
awk '{print $1}' "$3".tsv > motif.tsv
sort -n motif.tsv | uniq > motifs.sorted.tsv
sed 's/\//\./g' motifs.sorted.tsv > motifs.1.sorted.delete.tsv
sed 's/\"//g' motifs.1.sorted.delete.tsv > motifs.sorted.delete.tsv
dir=$(cd $(dirname $0);pwd)

grep "$3" $2 |grep "Qv" |awk -F'\t' '{print $9}'|sed 's/;coverage=/\tcoverage=/g' |awk '{print $2}'|sed 's/;frac=/\tfrac=/g'  |awk '{print $2";"$1}'|cut -f1,2,3 -d";" |sed 's/;coverage=/\t/g' |sed 's/frac=//g' |sed 's/;identificationQv=/\t/g' |awk '{print $3"\t"$1"\t"$2}'> "$3"_motif_basemods.txt

Rscript $dir/basemods.r "$3"_motif_basemods.txt
mv plot.pdf "$3"_plot.pdf 
mv fraction.plot.pdf "$3"_fraction.plot.pdf

####### step 1 methylation #######################################
for motif in `cat motifs.sorted.delete.tsv`
do
	grep ""$motif"\;cov" $2 > "$motif"_"$3"_1_motif
        grep "identificationQv" "$motif"_"$3"_1_motif > "$motif"_"$3"_motif
        cat "$motif"_"$3"_motif |cut -f9 |cut -f4 -d";" |sed -e 's/IPDRatio\=//g' > IPDR_"$motif".txt
        cat "$motif"_"$3"_motif|awk -F '\t' '{print $4}' > "$motif".bed
        paste "$motif".bed IPDR_"$motif".txt > IPDR_"$motif".bed.txt
        rm "$motif".bed
        rm IPDR_"$motif".txt
        less "$motif"_"$3"_motif |grep "Qv" |awk -F'\t' '{print $9}'|sed 's/;coverage=/\tcoverage=/g' |awk '{print $2}'|sed 's/;frac=/\tfrac=/g'  |awk '{print $2";"$1}'|cut -f1,2,3 -d";" |sed 's/;coverage=/\t/g' |sed 's/frac=//g' |sed 's/;identificationQv=/\t/g' |awk '{print $3"\t"$1"\t"$2}'> "$motif"_"$3"_motif_basemods.txt
        Rscript $dir/basemods.r "$motif"_"$3"_motif_basemods.txt
        mv plot.pdf "$motif"_"$3"_plot.pdf 
        mv fraction.plot.pdf "$motif"_"$3"_fraction.plot.pdf
done
cat IPDR*.bed.txt > whole.bed
rm *_1_motif

TESTIPDR=$(head whole.bed|awk '{print $2}'|grep "^[0-9]" )
if [ "$TESTIPDR" != "" ]
then
	echo"Wrong format of motifs.gff!"
        rm whole.bed
        rm *bed.txt
        exit 1
else

for motif in `cat motifs.sorted.delete.tsv`
do
        grep ""$motif"\;cov" $2 > "$motif"_"$3"_1_motif
        grep "identificationQv" "$motif"_"$3"_1_motif > "$motif"_"$3"_motif
        cat "$motif"_"$3"_motif | awk -F'\t' '{print $4"\t"$9}'|sed 's/;/\t/g' |awk '{print $1"\t"$9}' | sed -e 's/frac\=//g' > IPDR_"$motif".bed.txt
done
cat IPDR*.bed.txt > whole.bed
fi
rm *_1_motif
################methylation###################

for motif in `cat motifs.sorted.delete.tsv`
do
        less "$motif"_"$3"_motif|awk -F'\t' '{print $4"\t"$9}'|sed 's/;/\t/g'|awk '{print $1"\t"$5"\t"$10}'|sed 's/identificationQv=//g'|awk '{if($3>"'$5'"-1)print $0}'|sed 's/coverage=//g' |awk '{if($2> "'$6'"-1)print $1}'> "$motif"_"$3"_motif_methylation.gff.bed
done


#################step 2 undermethylation ####################################


for motif in `cat motifs.sorted.delete.tsv`
do
        less "$motif"_"$3"_motif|awk -F'\t' '{print $4"\t"$9}'|sed 's/;/\t/g'|awk '{print $1"\t"$5"\t"$9"\t"$10}'|sed 's/identificationQv=//g'|awk '{if($4>"'$5'"-1)print $0}'|sed 's/frac=//g' |awk '{if($3<"'$4'")print $0}' |sed 's/coverage=//g' |awk '{if($2>"'$6'"-1)print $1}'> "$motif"_"$3"_motif_undermethylation.gff.bed
done

#################step 3 no-methylation ####################################

for motif in `cat motifs.sorted.delete.tsv`
do
        echo $5
	less "$motif"_"$3"_motif|awk -F'\t' '{print $4"\t"$9}'|sed 's/;/\t/g'|awk '{print $1"\t"$5"\t"$9"\t"$10}'|sed 's/identificationQv=//g'|awk '{if ($4<='"$5"') print $0}'|sed 's/frac=//g' |sed 's/coverage=//g' |awk '{if($2>"'$6'"-1) print $1}'> "$motif"_"$3"_motif_unmethylation.gff.bed
	wc -l  "$motif"_"$3"_motif_unmethylation.gff*
done



rm motifs.sorted.delete.tsv
rm motifs.sorted.tsv
rm motif.tsv
rm test.tsv
rm motifs.1.sorted.delete.tsv 
rm "$3".tsv
rm *.bed.txt

