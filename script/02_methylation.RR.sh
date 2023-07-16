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
promoter=$4
cds=$5
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
# check methylation pattern

echo "head check methylation bed"
head motif.bed
echo "generate forward and reverse strand"
#####################1. prepare for regulation region##############################
grep '[[:space:]]+\+' genome.gff  > positive.gtf
grep '[[:space:]]-\+' genome.gff  > negative.gtf
dir=$(cd $(dirname $0);pwd)
Rscript $dir/get.tss.r $promoter $cds
echo "integrate regulaation region file"
awk '{print "'$sa'""\t"$1"\t"$2}' tss.po.bed > tss.bed
awk '{print "'$sa'""\t"$1"\t"$2}' tss.ne.bed >>tss.bed
echo "mapped in CDS region"
awk '{print "'$sa'""\t"$2-"'$cds'""\t"$2}' tss.po.bed > po.tss.cds.bed
awk '{print "'$sa'""\t"$1"\t"$1+"'$cds'"}' tss.ne.bed > ne.tss.cds.bed
echo "mapped in promoter region"
awk '{print "'$sa'""\t"$1"\t"$2-1-"'$cds'"}' tss.po.bed > po.tss.promoter.bed
awk '{print "'$sa'""\t"$1+1+"'$cds'""\t"$2}' tss.ne.bed >ne.tss.promoter.bed
awk '{if ($2<$3) print $0 }' po.tss.promoter.bed >po.promoter.bed
awk '{if ($2<$3) print $0 }' ne.tss.promoter.bed >ne.promoter.bed

#############2. select maped modification site in regulation region##################################

echo "do a bedtools merge on ${meth} and regulation region"
bedtools intersect -a tss.bed -b motif.bed -wa > tss.methylation.bed
bedtools intersect -a tss.bed -b motif.bed -wa -wb > tss.methylation.me.bed
wc -l tss.methylation.me.bed
bedtools intersect -a po.promoter.bed -b motif.bed -wa -wb> po.promoter.methylation.bed
bedtools intersect -a ne.promoter.bed -b motif.bed -wa -wb> ne.promoter.methylation.bed

countp1=`wc -l po.promoter.methylation.bed |awk '{print $1}'`

countp2=`wc -l ne.promoter.methylation.bed |awk '{print $1}'`
promotercount=$((${countp1}+${countp2}))
awk '{print $0"\t"$5-$3}'  po.promoter.methylation.bed > po.promoter.methylation.txt
awk '{print $0"\t"$2-$5}'  ne.promoter.methylation.bed > ne.promoter.methylation.txt


bedtools intersect -a po.tss.cds.bed -b motif.bed -wa -wb > po.cds.methylation.bed
bedtools intersect -a ne.tss.cds.bed -b motif.bed -wa -wb > ne.cds.methylation.bed
###################################CDS##################################

wc -l po.cds.methylation.bed
wc -l ne.cds.methylation.bed
count1=`wc -l po.cds.methylation.bed |awk '{print $1}'`
count2=` wc -l ne.cds.methylation.bed |awk '{print $1}'`
count=$((${count1}+${count2}))
awk '{print $0"\t"$5-$2}' po.cds.methylation.bed > po.cds.methylation.txt
awk '{print $0"\t"$3-$5}' ne.cds.methylation.bed > ne.cds.methylation.txt

echo "generate a fatsa file means this regulation methylated"
#sort -n tss.methylation.bed | uniq > tss.methylation.sorted.bed
#integrate cds and regulation bed file
#cat tss.methylation.sorted.bed coding.methylation.sorted.bed > tss.methylation.sorted.1.bed
bedtools getfasta -fi genome.fna  -bed tss.methylation.bed -fo $SAMPLE$meth.fimo.fasta
cp $SAMPLE$meth.fimo.fasta "$1"_methylationGene.fasta
echo "find methylated genes"


awk '{print $5}' negative.gtf > test.bed
awk '{print $4}' positive.gtf >> test.bed
##############################3. output the gene decription##########################################


cp po.cds.methylation.txt ./cds.txt
awk '{print $2}' cds.txt > cds.tss.bed
awk '{print $3"\t"$4}' positive.gtf > test.bed
grep -wf cds.tss.bed test.bed > over.cds.txt
grep -wf over.cds.txt positive.gtf >gene.po.cds.gtf
cat gene.po.cds.gtf |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp.cds
cp locus.tmp.cds "$1"_rr_locustag.txt
grep -wf locus.tmp.cds genome.tsv |cut -f1,4,7 > gene.cds.list.txt
awk '{print $4"\t"$5"\t"$7}' gene.po.cds.gtf > tmp.txt
paste tmp.txt gene.cds.list.txt > gene.txt
awk '{print $2"\t"$3"\t"$5"\t"$7}' po.cds.methylation.txt >po.cds.methylation.1.txt
awk '{print $1"\t"$0}' po.cds.methylation.1.txt >test.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $2"\t"$3"\t"$4"\t"$5"\t"a[$2]}' gene.txt test.txt > po.cds.me.txt
cp locus.tmp.cds cds.locus.txt

cp ne.cds.methylation.txt ./cds.txt 
awk '{print $3}' cds.txt > cds.tss.bed
awk '{print $5"\t"$6}' negative.gtf > test.bed
grep -wf cds.tss.bed test.bed > over.cds.txt
grep -wf over.cds.txt negative.gtf > gene.ne.cds.gtf
cat gene.ne.cds.gtf|cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp.cds

cat locus.tmp.cds >> "$1"_rr_locustag.txt

grep -wf locus.tmp.cds genome.tsv |cut -f1,4,7 > gene.cds.list.txt
cat cds.locus.txt locus.tmp.cds > cds.locus_1.txt
mv cds.locus_1.txt cds.locus.txt
awk '{print $4"\t"$5"\t"$7}' gene.ne.cds.gtf >tmp.txt
paste tmp.txt  gene.cds.list.txt > gene.txt
awk '{print $3"\t"$5"\t"$7}' ne.cds.methylation.txt > ne.cds.methylation.1.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$2]"\t"$0}'  ne.cds.methylation.1.txt gene.txt|awk '{print $1-"'$cds'""\t"$0}' > ne.cds.me.txt
cat po.cds.me.txt ne.cds.me.txt|awk '{print "'$SAMPLE'""\t""'$meth'""\t""CDS""\t"$0}'> cds.gene.me.txt

###################################regulation##################################


cp po.promoter.methylation.txt ./promoter.methylation.bed
awk '{print $3+1}' promoter.methylation.bed > promoter.tss.bed
awk '{print $3"\t"$4}' positive.gtf > test.bed
grep -wf promoter.tss.bed  test.bed > over.promoter.txt
grep -wf over.promoter.txt positive.gtf >gene.po.promoter.gtf
cat gene.po.promoter.gtf |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp
cat locus.tmp.cds >> "$1"_rr_locustag.txt
grep -wf locus.tmp genome.tsv |cut -f1,4,7 > gene.promoter.list.txt
awk '{print $4"\t"$5"\t"$7}' gene.po.promoter.gtf >tmp.txt
paste tmp.txt gene.promoter.list.txt > gene.txt
awk '{print $2"\t"$3+1"\t"$5"\t"$7}' po.promoter.methylation.txt >  po.promoter.methylation.1.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $1"\t"$2"\t"$3"\t"$4"\t"a[$2]}' gene.txt po.promoter.methylation.1.txt > po.promoter.me.txt
cp locus.tmp promoter_locus.txt
cp ne.promoter.methylation.txt ./promoter.methylation.bed
awk '{print $2-1}' promoter.methylation.bed > promoter.tss.bed
awk '{print $5"\t"$6}' negative.gtf > test.bed
grep -wf promoter.tss.bed  test.bed > over.promoter.txt
grep -wf over.promoter.txt negative.gtf >gene.ne.promoter.gtf
cat gene.ne.promoter.gtf |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp
cat locus.tmp.cds >> "$1"_rr_locustag.txt
grep -wf locus.tmp genome.tsv |cut -f1,4,7 > gene.promoter.list.txt
awk '{print $4"\t"$5"\t"$7}' gene.ne.promoter.gtf >tmp.txt
paste tmp.txt gene.promoter.list.txt > gene.txt
awk '{print $2-1"\t"$3"\t"$5"\t"$7}' ne.promoter.methylation.txt >  ne.promoter.methylation.1.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$2]"\t"$0}' ne.promoter.methylation.1.txt gene.txt > ne.promoter.me.txt
cat promoter_locus.txt locus.tmp > promoter_locus_1.txt
mv promoter_locus_1.txt promoter_locus.txt
cat po.promoter.me.txt ne.promoter.me.txt|awk '{print"'$SAMPLE'""\t""'$meth'""\t""promoter""\t"$0}' > promoter.gene.me.txt



#statistic
echo "Strain    Methylation     nCDS    nPROMOTER       nRR" > "$1"_result.txt
echo "$SAMPLE   $meth   $count       $promotercount  $(($count+$promotercount))" >> "$1"_result.txt
cat cds.gene.me.txt promoter.gene.me.txt > "$1"_methylationGene.txt
cp "$1"_methylationGene.txt "$1"_raw_methylationGene.txt

awk '{print $6"\t"$0}' "$1"_methylationGene.txt > test.txt
awk '{print $2"\t"$1}' whole.bed > ipdr.mid.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$2]}'  test.txt ipdr.mid.txt > methylationGene.txt
grep "$SAMPLE" methylationGene.txt >methylationGene.1.txt

echo "Fraction	 	 	Strain	Methylation	Region	RRS	RRE	Methsite	distance	start	end	strand	Locus_tag    gene name	description"> "$1"_methylationGene.1.txt
grep "$SAMPLE" methylationGene.txt >> "$1"_methylationGene.1.txt

cat "$1"_methylationGene.1.txt|awk -v FS="\t" -v OFS="\t" '{print $4,$5,$6,$7,$8,$9,$1,$10,$11,$12,$13,$14,$15,$16}' > "$1"_methylationGene.txt

rm methylationGene.txt
rm ipdr.mid.txt
rm "$1"_methylationGene.1.txt
rm tss.bed
rm tss.po.bed
rm tss.ne.bed
rm tmp.gff
rm cds.gene.me.txt
rm cds.tss.bed
rm cds.txt
rm gene.cds.list.txt
rm gene.ne.cds.gtf
rm gene.ne.promoter.gtf
rm gene.po.cds.gtf
rm gene.po.promoter.gtf
rm gene.promoter.list.txt
rm gene.txt
rm genome.fna
rm genome.fna.fai
rm genome.gff
rm genome.tsv
rm locus.tmp
rm locus.tmp.cds
rm motif.bed
rm ne.cds.methylation.1.txt
rm ne.cds.methylation.bed
rm ne.cds.methylation.txt
rm ne.cds.me.txt
rm negative.gtf
rm ne.promoter.bed
rm ne.promoter.methylation.1.txt
rm ne.promoter.methylation.bed
rm ne.promoter.methylation.txt
rm ne.promoter.me.txt
rm ne.tss.cds.bed
rm ne.tss.promoter.bed
rm over.cds.txt
rm over.promoter.txt
rm po.cds.methylation.1.txt
rm po.cds.methylation.bed
rm po.cds.methylation.txt
rm po.cds.me.txt
rm po.promoter.bed
rm po.promoter.methylation.1.txt
rm po.promoter.methylation.bed
rm po.promoter.methylation.txt
rm po.promoter.me.txt
rm positive.gtf
rm po.tss.cds.bed
rm po.tss.promoter.bed
rm promoter.gene.me.txt
rm promoter.methylation.bed
rm promoter.tss.bed
rm test.bed
rm test.txt
rm tmp.txt
rm tss.methylation.bed
rm tss.methylation.me.bed
#rm tss.methylation.sorted.bed
rm "$1"_raw_methylationGene.txt
#rm "$1"_methylationGene.1.txt








