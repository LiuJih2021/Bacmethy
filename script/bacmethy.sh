#!/bin/bash -login
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -j y
# Time-stamp: <2021-06-27 16:49:42 LiuJihong>
usage() {
    echo "
Bacmethy V1.0.1 01/01/2022
<liujihong2021@126.com>
usage:  bacmethy.sh  -m <motifs.gff> -s <motifs.csv> -g <genome.fa> -p <prefix> -t <m6A, m4C or m5C>  [options] -T  -d -i -c -b -a -r 
requirement: locally installed PROKKA,  bedtools and MEME  softwares


options
    -m    FILE    a motif detected file in gff format from PacBio portal (required)
    -s    FILE    a motifs summary file (required)
    -g      FILE    a genome file which complete sequenced (required)
    -p      STRING  prefix of output (required; usually a strains name or a sample name)
    -t      STRING  a type of methylation type (required; one of m6A, m4C or m5C)
    -T      FOLDER  a floder contains TFs files in meme format (required; users can get from Calcute_PPM_console_final.py)
    -d      FLOAT   a undermethylation thresholds of fraction (default: 0.75)
    -i      INT     a identificationQV thresholds (default: 40)
    -c      INT     a coverage thresholds (default: 30)
    -b      INT     number of bps before TSS (default: 500)
    -a      INT     number of bps after TSS (default: 100)
    -r      FILE    FASTA or GBK file to use as 1st priority (default '')
    "
}

while getopts ":hm:s:g:p:t:T:d:b:a:r:i:n:c:" opt
do
    case $opt in
        h)  usage
            exit 0
            ;;
        m)    METHYLATIONGFF=$OPTARG ;;
        s)    METHYLATIONCSV=$OPTARG ;;
        g)  FASTA=$OPTARG ;;
        p)  sample=$OPTARG ;;
        t)  METHYLATIONTYPE=$OPTARG ;;
        T) MEME=$OPTARG ;;
        d)  undermethylation_threshold=$OPTARG ;;
        b)  promoter=$OPTARG ;;
        a)  cds=$OPTARG ;;
        r)  proteins=$OPTARG ;;
        i)  identification=$OPTARG ;;
        n)  npu=$OPTARG ;;
        c) coverage=$OPTARG ;;
        :)
        echo "Option -$OPTARG requires an argument."
	    usage
        exit 1
        ;;
	    ?) 
	        echo "Invalid option: -$OPTARG"
	        usage      
            exit 1
	    ;;
    esac    
done

#If options of 'd', 'c'and 'a' are not defined, setting to default values
undermethylation_threshold=${undermethylation_threshold:-0.75}
promoter=${promoter:-500}
cds=${cds:-100}
identification=${identification:-40}
coverage=${coverage:-30}
npu=${npu:-5}



#check locally software
PROKKA=`which prokka || true`
if [[ -z "$PROKKA" ]]
then
    usage
	echo "PROKKA not found"
	exit 1
fi
bedtools=`which bedtools || true`
if [[ -z "$bedtools" ]]
then
    usage
	echo "bedtools not found"
	exit 1
fi
fimo=`which fimo || true`
if [[ -z "$fimo" ]]
then
    usage
	echo "fimo not found"
	exit 1
fi
####### running btpmi ########
start_time=`date +%s`
dir=$(cd $(dirname $0);pwd)
mkdir -p $dir/tmp
tmp_fifofile="$dir/tmp/$$.fifo"

mkfifo $tmp_fifofile
exec 6<>$tmp_fifofile
rm $tmp_fifofile


for((i=0;i<${npu};i++));do
    echo
done >&6


if [[ ! -d "./methylation/" ]]
then
    mkdir methylation
    mkdir undermethylation
    mkdir unmethylation
fi


    echo "variables to process"
    echo $sample $species $METHYLATIONGFF $METHYLATIONTYPE  $FASTA $promoter $cds $undermethylation_threshold ${npu} 
    echo $coverage
    echo "running prokka,"
    if [[ -r "$proteins" ]]
    then
    prokka $FASTA --outdir ${sample}  --prefix ${sample}  --kingdom Bacteria  --proteins $proteins 
    else
    prokka $FASTA --outdir ${sample}  --prefix ${sample}  --kingdom Bacteria
    fi
    
    echo "running different methylation type"
    echo "The undermethylation threshold is "$undermethylation_threshold""
    
    $dir/01_methylation_filter.sh $METHYLATIONCSV $METHYLATIONGFF $METHYLATIONTYPE $undermethylation_threshold $identification $coverage
####################################QC####################################



motifQ=`wc -l *motif|tail -n -1|awk '{print $1}'`
methQ=`wc -l *methylation.gff.bed|tail -n -1|awk '{print $1}'`
QC=`echo "scale=3; $methQ / $motifQ "|bc`
if [ $(echo "$QC < 0.5" | bc) -eq 1 ]
then
    echo "low quality!!!"

fi


###################step 2 ###########################

motifGff=`ls -1 *motif`
for motifgff in $motifGff
do
    #############methylated#######################
    echo "running methylated genes"
    awk '{print "'"$METHYLATIONTYPE"'\t"$1}' "$motifgff"_methylation.gff.bed > "$motifgff"_methylation.gff.bed.1   
    grep -wf "$motifgff"_methylation.gff.bed.1 "$motifgff" > "$motifgff".methylation

#    grep -wf "$motifgff"_methylation.gff.bed "$motifgff" > "$motifgff".methylation
    $dir/02_methylation.RR.sh  $motifgff.methylation  ${sample} $METHYLATIONTYPE $promoter $cds
    
    if [ -n "$MEME" ]; then
    meme_file=`ls -1 $MEME`
    for memefile in $meme_file
    do
        echo "running FIMO scan"
        fimo $MEME/$memefile  ${sample}$METHYLATIONTYPE.fimo.fasta
        echo "integrated genes"
        $dir/03_afterfimo.sh ./fimo_out/*.gff ${sample} $motifgff.methylation $memefile $promoter $cds
        rm -r fimo_out
    done
    fi
    mkdir -p ./methylation/"$motifgff"_RR
    mv "$motifgff".methylation* ./methylation/"$motifgff"_RR

    cp ./methylation/"$motifgff"_RR/"$motifgff".methylation ./
    $dir/02_methylation.gene.sh  $motifgff.methylation  ${sample} $METHYLATIONTYPE $promoter $cds
    mkdir -p ./methylation/"$motifgff"_CDS
    mv "$motifgff".methylation* ./methylation/"$motifgff"_CDS


    #############undermethylated#######################
    echo "running undermethylated genes"
    
    awk '{print "'"$METHYLATIONTYPE"'\t"$1}' "$motifgff"_undermethylation.gff.bed > "$motifgff"_undermethylation.gff.bed.1
    grep -wf "$motifgff"_undermethylation.gff.bed.1 "$motifgff" > "$motifgff".undermethylation

    #grep -wf "$motifgff"_undermethylation.gff.bed "$motifgff" > "$motifgff".undermethylation
    
    $dir/02_methylation.RR.sh  "$motifgff".undermethylation  ${sample} $METHYLATIONTYPE $promoter $cds
    
    if [ -n "$MEME" ]; then
    meme_file=`ls -1 $MEME`
    for memefile in $meme_file
    do
        echo "running FIMO scan"
        fimo $MEME/$memefile  ${sample}$METHYLATIONTYPE.fimo.fasta

        echo "integrated genes"
         $dir/03_afterfimo.sh ./fimo_out/*.gff ${sample} $motifgff.undermethylation $memefile $promoter $cds
         rm -r fimo_out
    done
    fi
    mkdir -p ./undermethylation/"$motifgff"undermethylation_RR
    mv "$motifgff".undermethylation* ./undermethylation/"$motifgff"undermethylation_RR

    cp ./undermethylation/"$motifgff"undermethylation_RR/"$motifgff".undermethylation ./
    $dir/02_methylation.gene.sh  "$motifgff".undermethylation  ${sample} $METHYLATIONTYPE $promoter $cds
    
    mkdir -p ./undermethylation/"$motifgff"undermethylation_CDS
    mv "$motifgff".undermethylation* ./undermethylation/"$motifgff"undermethylation_CDS


   
    #############unmethylated#######################
    echo "running unmethylated genes"
    awk '{print "'"$METHYLATIONTYPE"'\t"$1}' "$motifgff"_unmethylation.gff.bed > "$motifgff"_unmethylation.gff.bed.1 
    grep -wf "$motifgff"_unmethylation.gff.bed.1 "$motifgff" > "$motifgff".unmethylation
#    grep -wf "$motifgff"_unmethylation.gff.bed "$motifgff" > "$motifgff".unmethylation
    $dir/02_methylation.RR.sh  "$motifgff".unmethylation  ${sample} $METHYLATIONTYPE $promoter $cds
    
    if [ -n "$MEME" ]; then
    meme_file=`ls -1 $MEME`
    for memefile in $meme_file
    do
        echo "running FIMO scan"
        fimo $MEME/$memefile  ${sample}$METHYLATIONTYPE.fimo.fasta

        echo "integrated genes"
        $dir/03_afterfimo.sh ./fimo_out/*.gff ${sample} $motifgff.unmethylation $memefile $promoter $cds
         rm -r fimo_out
    done
    fi
    mkdir -p ./unmethylation/"$motifgff"unmethylation_RR
    mv "$motifgff".unmethylation* ./unmethylation/"$motifgff"unmethylation_RR

    cp ./unmethylation/"$motifgff"unmethylation_RR/"$motifgff".unmethylation ./
    $dir/02_methylation.gene.sh  "$motifgff".unmethylation  ${sample} $METHYLATIONTYPE $promoter $cds
    
    
    mkdir -p ./unmethylation/"$motifgff"unmethylation_CDS
    mv "$motifgff".unmethylation* ./unmethylation/"$motifgff"unmethylation_CDS
    read -u6
    {
        echo "$motifgff finished!"
        echo >&6
    } &
done 

wait

rm whole.bed
rm cds.locus.txt
rm *fimo.fasta
rm promoter_locus.txt

stop_time=`date +%s`
echo "Time:`expr $stop_time - $start_time`"

exec 6>&-
echo "over !"

rm -r $dir/tmp
rm *gff.bed*
