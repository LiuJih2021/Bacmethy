#!/usr/bin/env bash
#bedtools v2.30.0
# Time-stamp: <2021年 07月 01日 星期四 11:38:48 CST liujihong>

usage() {
    echo "usage: bash 02_methylation.gene.sh <METHYLATIONGFF> <SAMPLE> "
    echo "where: <METHYLATIONGFF> is a specific DNA methylation file"
    echo "       <SAMPLE> is the strain&DNA methylation name"
    echo "This script does the job of fetching the list of variables to pass to"
    echo "the script a1_BTPMI.sh"
}

# Minimal argument checking

if [ $# -lt 1 ]; then
    usage
    exit
fi
# Set variable for input file
SAMPLE=$2
meth=$3
echo "Start"
date
#get methylatioin site bed
cp $1 tmp.gff
cp $SAMPLE/*.gff genome.gff
cp $SAMPLE/*.fna genome.fna
cp $SAMPLE/*.tsv genome.tsv
sa=`grep ">" genome.fna |sed 's/^.//'`

echo "check genome name $sa"
#sed -i '1,3d' tmp.gff
cat tmp.gff |awk -F '\t' '{print "'$sa'""\t"$4"\t"$4}' > motif.bed

echo "head check methylation bed"
head motif.bed
echo "generate forward and reverse strand"
#positive bed
grep '[[:space:]]+\+' genome.gff  > positive.gtf
grep '[[:space:]]-\+' genome.gff  > negative.gtf

################################coding region################################

echo "mapped in CDS region"
awk '{print "'$sa'""\t"$4"\t"$5}' positive.gtf >codingRegion.bed

bedtools intersect -a codingRegion.bed -b motif.bed -wa > coding.methylation.bed
bedtools intersect -a codingRegion.bed -b motif.bed -wa -wb >coding.methylation.me.bed
#wc -l coding.methylation.me.bed

awk '{print $0"\t"$5-$2}'  coding.methylation.me.bed  > coding.methylation.1.me.txt

awk '{print "'$sa'""\t"$4"\t"$5}' negative.gtf >codingRegion.bed
bedtools intersect -a codingRegion.bed -b motif.bed -wa >> coding.methylation.bed
bedtools intersect -a codingRegion.bed -b motif.bed -wa -wb >coding.methylation.2.me.bed
awk '{print $0"\t"$3-$5}'  coding.methylation.2.me.bed  > coding.methylation.2.me.txt

cat coding.methylation.1.me.txt coding.methylation.2.me.txt >coding.methylation.me.txt 



sort -n coding.methylation.bed |uniq >coding.methylation.sorted.bed
#bedtools getfasta -fi genome.fna  -bed coding.methylation.sorted.bed  -fo  $SAMPLE$meth.fimo.fasta
wc -l  coding.methylation.sorted.bed


##############################coding region methylation gene##########################################
cp coding.methylation.me.txt ./cds.txt

awk '{print $2"\t"$3}' cds.txt > cds.tss.bed
grep -wf cds.tss.bed  genome.gff >gene.po.cds.gtf
cat gene.po.cds.gtf |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp.cds
grep -wf locus.tmp.cds genome.tsv |cut -f1,4,7 > gene.cds.list.txt
awk '{print $4"\t"$5"\t"$7}' gene.po.cds.gtf > tmp.txt
paste tmp.txt gene.cds.list.txt > gene.txt
awk '{print $2"\t"$3"\t"$5"\t"$7}' coding.methylation.me.txt >po.cds.methylation.1.txt
awk '{print $1"\t"$0}' po.cds.methylation.1.txt >test.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $2"\t"$3"\t"$4"\t"$5"\t"a[$2]}' gene.txt test.txt > coding.cds.me.txt
cat coding.cds.me.txt |awk '{print "'$SAMPLE'""\t""'$meth'""\t""CodingRegion""\t"$0}'> coding.gene.me.txt



#statistic
#cdscount=`wc -l cds.gene.me.txt|awk '{print $1}'`
#promotercount=`wc -l promoter.gene.me.txt|awk '{print $1}'`
Codingcount=`wc -l coding.methylation.sorted.bed |awk '{print $1}'`
echo "Strain    Methylation nCodingRegionCount" > "$1"_result.txt
echo "$SAMPLE   $meth   $Codingcount" >> "$1"_result.txt
#awk '{print "'$SAMPLE'""\t""'$meth'""\t""'$cdscount'""\t""'$promotercount'"}' >>result.txt
cp coding.gene.me.txt "$1"_methylationGene.txt


awk '{print $6"\t"$0}' "$1"_methylationGene.txt > test.txt
awk '{print $2"\t"$1}' whole.bed > ipdr.mid.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$2]}'  test.txt ipdr.mid.txt > methylationGene.txt
grep "$SAMPLE" methylationGene.txt >methylationGene.1.txt

echo "Fraction	 	 	Strain	Methylation	Region	RRS	RRE	Methsite	distance	start	end	strand	Locud_tag    gene name	description"> "$1"_methylationGene.1.txt
grep "$SAMPLE" methylationGene.txt >> "$1"_methylationGene.1.txt

cat "$1"_methylationGene.1.txt|awk -v FS="\t" -v OFS="\t" '{print $4,$5,$6,$7,$8,$9,$1,$10,$11,$12,$13,$14,$15,$16}' > "$1"_methylationGene.txt

rm cds.tss.bed
rm cds.txt
rm coding.cds.me.txt
rm coding.gene.me.txt
rm coding.methylation.1.me.txt
rm coding.methylation.2.me.bed
rm coding.methylation.2.me.txt
rm coding.methylation.bed
rm coding.methylation.me.bed
rm coding.methylation.me.txt
rm coding.methylation.sorted.bed
rm codingRegion.bed
rm gene.cds.list.txt
rm gene.po.cds.gtf
rm gene.txt
rm genome.fna
rm genome.gff
rm genome.tsv
rm ipdr.mid.txt
rm locus.tmp.cds
rm methylationGene.1.txt
rm methylationGene.txt
rm motif.bed
rm negative.gtf
rm po.cds.methylation.1.txt
rm positive.gtf
rm test.txt
rm tmp.gff
rm tmp.txt
rm "$1"_methylationGene.1.txt





