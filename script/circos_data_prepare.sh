#usage bash /PATH/TO/Circos_data_prepare.sh PREFIX /PATH/TO/SCRIPT/FOLDER/
mkdir circos_tmp
cp $2*.conf ./circos_tmp
cd circos_tmp
SAMPLE=$1
GENOME_LENGTH=$(seqkit stat ../$SAMPLE/$SAMPLE.fna |awk '{print $5}'|sed -e 1d)
echo chr - chr1 chr1 0 $GENOME_LENGTH set2-5-qual-3 |sed -e 's/,//g'> karyotype.txt
less ../$SAMPLE/$SAMPLE.gff |grep '[[:space:]]+\+'|awk -F"\t" '{print "chr1""\t"$4"\t"$5"\t"10}'|sort -n |uniq   > $SAMPLE.positive.bed
less ../$SAMPLE/$SAMPLE.gff |grep '[[:space:]]-\+'|awk -F"\t" '{print "chr1""\t"$4"\t"$5"\t"10}'|sort -n |uniq   > $SAMPLE.negative.bed
cat ../methylation/*RR/*_motif.methylation |grep '[[:space:]]+\+'|awk -F"\t" '{print "chr1""\t"$4"\t"$5+5"\t"10}' |sort -n |uniq > circos_4_value.txt
cat ../methylation/*RR/*_motif.methylation |grep '[[:space:]]-\+'|awk -F"\t" '{print "chr1""\t"$4"\t"$5+5"\t"10}'|sort -n |uniq  > circos_6_value.txt
cat ../un*thylation/*RR/*un*tion |grep '[[:space:]]+\+'|awk -F"\t" '{print "chr1""\t"$4"\t"$5"\t"10}'|sort -n |uniq  > circos_51_value.txt
cat ../un*thylation/*RR/*un*tion |grep '[[:space:]]-\+'|awk -F"\t" '{print "chr1""\t"$4"\t"$5"\t"10}'|sort -n |uniq  > circos_52_value.txt
cat circos_51_value.txt circos_52_value.txt > circos_5_value.txt
cat circos_4_value.txt circos_6_value.txt > circos_7_value.txt
cat $SAMPLE.positive.bed |awk '{print $0"\t""fill_color=chr"int(16*rand())}' |sort -n |uniq > circos_2_highlight.bed
cat $SAMPLE.negative.bed |awk '{print $0"\t""fill_color=chr"int(16*rand())}' |sort -n |uniq > circos_3_highlight.bed
circos -conf circos_set.conf
cd ..
