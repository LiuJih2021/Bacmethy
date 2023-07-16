#!/usr/bin/env bash
SAMPLE=$2
promoter=$3
cds=$4
cp $1 ./fimo.gff
cp $SAMPLE/*.gff genome.gff
cp $SAMPLE/*.fna genome.fna
cp $SAMPLE/*.tsv genome.tsv
sa=`grep ">" genome.fna |sed 's/^.//'`
awk '{print $1}' fimo.gff |sed -e "s/$sa://g"|sed -e 's/-/\t/g' >tmp.bed
#awk '{print $4"\t"$5"\t"$6}' fimo.gff > score.tmp.txt
awk '{print $4"\t"$5"\t"$10}' fimo.gff |cut -f1 -d";" > score.tmp.txt
paste tmp.bed score.tmp.txt > tmp.bed.txt
sed -i '1d' tmp.bed
awk '{print $1+"'$cds'"}' tmp.bed >tss.bed
awk '{print $2-"'$cds'"}' tmp.bed >>tss.bed
grep -wf tss.bed genome.gff  |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp
grep -wf locus.tmp genome.tsv |cut -f4,7 >gene.list.txt
awk '{print $2"\t"$0}' tmp.bed.txt >tmp.1.txt
awk '{print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6}' tmp.1.txt >tmp.2.txt
cat tmp.2.txt >>tmp.1.txt
awk '{print $4"\t"$0}' $3_methylationGene.txt > gene.txt
awk '{print $5"\t"$0}' $3_methylationGene.txt >> gene.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$2]}' gene.txt tmp.1.txt  > tmp.txt
echo "Strain	Methylation	TF binding start	TF binding end	FIMO qvalue	Region	RRS	RRE	Methsite	Fraction	distance	start	end	strand  Locus_tag	gene name	description"> "$3"_"$4"_binding_methylation.txt
grep "$SAMPLE" tmp.txt > "$3"_TF.txt
cat "$3"_TF.txt|awk -v FS="\t" -v OFS="\t" '{$4=$4+$2;$5=$5+$2;$1="";$2="";$3="";$7="";print $8,$9,$4,$5,$6,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' >> "$3"_"$4"_binding_methylation.tmp.txt
sort -n "$3"_"$4"_binding_methylation.tmp.txt |uniq >> "$3"_"$4"_binding_methylation.txt
rm "$3"_"$4"_binding_methylation.tmp.txt
#sed 's/\_.//' gene.list.txt| sort -n|uniq >gene.1.list.txt
rm tmp.bed
rm tss.bed
rm locus.tmp
rm tmp.1.txt
rm tmp.bed.txt
rm tmp.txt
rm score.tmp.txt
rm gene.txt
rm gene.list.txt
rm tmp.2.txt
rm fimo.gff	
