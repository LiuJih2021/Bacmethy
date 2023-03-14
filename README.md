Bacmethy V 1.0.1
=================
Summary
------------------

Bacmethy, a customizable pipeline based on Docker, reliably calculates the enrichment significance of methylated and un(der)methylated motifs in regulatory and coding regions of the genome, identifies the DNA methylation and transcriptional factor binding co-effected genes, facilitating the prediction of regulation effects on DNA methylation, and achieving a more comprehensive understanding of bacterial epigenomes.<br>


Requirment
---------------

1. [PROKKA](https://github.com/tseemann/prokka "PROKKA Github")<br>
2. [bedtools](https://bedtools.readthedocs.io/en/latest/ "bedtools Documentation")<br>
3. [meme](https://github.com/cinquin/MEME "meme Github")<br>


Installation
---------------
### 1. WEB Sever
[Bacmethy Website](https://Bacmethy.med.sustech.edu.cn)
The website includes methylation analysis and TFs binding prediction of Bacmethy. <br>
We recommend users visit [Bacmethy](https://Bacmethy.med.sustech.edu.cn) to get result easily.<br>
 
### 2. docker users
#### for Mac or Windows users
With docker developing, it become more easy to use Linux software in other operation system, like Mac or Windows. Users can install [**docker Desktop** ](https://www.docker.com/get-started/) and use Bacmethy in Mac or Windows. <br>
#### First Download
1. After [docker](https://docs.docker.com/engine/install/) installation. <br>
    docker pull liujihong/Bacmethy:2.0
    docker run  -t -i liujihong/Bacmethy:2.0 /bin/bash
2. re-enter docker. <br>
    docker start <container_ID>
    docker attach <container_ID>

### 3. others
1. Download the script to your working directory. <br>
2. install required softwares through Conda. <br>
#### First install conda 
    wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    bash Miniconda2-latest-Linux-x86_64.sh
#### creat a Bacmethy conda environment
    conda -n Bacmethy
    conda activate Bacmethy
#### conda install required softwares
    conda install -c bioconda prokka
    conda install -c bioconda bedtools
    conda install -c bioconda meme

Run Bacmethy
-----------------
### 1. Prepare motifs file
Use Base Modification Analysis Application in [SMRT-LINK](https://www.pacb.com/support/software-downloads/) to identify common bacterial base modifications, and then optionally analyze the methyltransferase recognition motifs. Detection can use an in-silico control consisting of expected kinetic signals.<br>
The result of 'Base Modification Analysis ' (motifs.gff and motif.csv)are required in Bacmethy analysis. <br>
1. motifs.gff: Compliant with the GFF Version3 [specification](http://www.sequenceontology.org/gff3.shtml). Each template position/strand pair whose probability value exceeds the probability value threshold appears as a row. The template position is 1-based, per the GFF specifications. The strand column refers to the strand carrying the detected modification, which is the opposite strand from those used to detect the modification. The GFF confidence column is a Phred-transformed probability value of detection.<br>
2. motifs.csv: Containsonerowforeach(referenceposition, strand) pair that appeared in the Data Set with coverage at least x.
x defaults to 3, but is configurable with the --minCoverage option. The reference position index is 1-based for compatibility with the GFF file in the R environment. Note that this output type scales poorly and is not recommended for large genomes; the HDF5 output should perform much better in these cases.
ref: https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v10.1.pdf
### 2. Prepare genome file 
PacBio provide HGAP4(Hierarchical Genome Assembly Process) in [SMRT-LINK](https://www.pacb.com/support/software-downloads/) to generate high quality de novo assemblies of genomes, using Continuous Long Reads.<br>

**Complete genome** are required which **from same strains** in step 1 Prepare motifs files. 
### 3. Prepare meme file [OPTION]
TFs binding predictiion as an optional function module of Bacmethy. The TF matrix file must in meme format.
### 4. Run Bacmethy.sh
    bash Bacmethy.sh -m motifs.gff -s motifs.csv -g genome.fasta -p PREFIX -t m6A 
#### Usage
    bash Bacmethy.sh -m <motifs.gff> -s <motifs.csv> -g <genome.fa> -p <prefix> -t <m6A, m4C or m5C> [options] -d -b -a -r -T -n
    requirement: locally installed PROKKA,  bedtools and MEME  softwares


    required
        -m      FILE    a motif detected file in gff format from PacBio portal (required)
        -s      FILE    a motifs summary file (required)
        -g      FILE    a genome file which complete sequenced (required)
        -p      STRING  prefix of output (required; usually a strains name or a sample name)
        -t      STRING  a type of methylation type (required; one of m6A, m4C or m5C)
       
    options
        -T      FOLDER  a floder contains TFs files in meme format (required; users can get from Calcute_PPM_console_final.py)
        -d      FLOAT   a undermethylation thresholds of fraction (default: 0.75)
        -i      FLOAT   a unmethylation thresholds of identificationQv (default: 40)
        -c      FLOAT   a unmethylation thresholds of coverage (default: 30)
        -b      INT     number of bps before TSS (default: 500)
        -a      INT     number of bps after TSS (default: 100)
        -r      FILE    FASTA or GBK file to use as 1st priority (default '')
        -n      INT     Number of CPUs to use (default '8')
#### Notice
The genome must be the same file which used methylation motif detection analysis.<br>
#### option: **-r** add standard strains reference
To better annotate bacterial genome, we recommend users use **-r proteins.faa** which from a standard strain in same species. <br>
#### option: **-d** methylation level direction
The epigenetic regulation involved in DNA methylation is complex and dynamic. All kinds of regulation on DNA may affect DNA methylation events, resulting in different signals detected by sequencing machine, like PacBio, Nanopore... and classified different levels of DNA methylation events. In order to ensure the authenticity of methylation events, we require the sequencing coverage to be higher than **30%**. According to different signals of DNA methylation, the software is divided into three methylation levels (methylation, undermethylation and unmethylation) for analysis respectively, of which the default thresholds are:<br>
1. methylation: identificationQV > 40 <br>
2. undermethylation: identificationQV > 40 && fraction < 0.75 <br>
3. unmethylation: identifictionQV < 40 <br>
Bacmethy set a default threshold of undermethylation(fraction < 0.75) and identificationQV > 40, users also can set user-defined threshold by using parameters **-d** to change  fraction and **-i** to change identificationQV. <br>
### 5. example
    Bacmethy.sh -m motifs.gff -s motif_summary.csv -g genomic.fna -p K12 -t m6A -T /Your/Path/To/TF/meme

### 6. Running Time
Bacmethy uses parallel processing to decrease running time on multicore computers.  users can set Running GPU by parameter **-n**. 

### 7. Output
structure of output files
-   methylation
    -   motif_CDS
    -   motif_RR
-   undermethylation
    -   motif_CDS
    -   motif_RR
-   unmethylation
    -   motif_CDS
    -   motif_RR
#### Gene Features
There are 3 gene features, promoter(default: 500bp before TSS), CDS_RR(default: 100bp after TSS), Coding Region(Whole Gene coding region). And we defined both promoter and CDS as Regulation Region(RR) in Bacmethy. <br>
##### difference between CDS_RR and Coding Region<br>
Gene transcription initiation is a multi regulated process. So the start of whole gene body also may play a role in gene transcription initiation regulation, then the **CDS_RR** was used to descript that methylation sites happened at the start CDS region which regulated gene transcription. <br>
Somes users may focus on the methylation sites happened at gene body. Bacmethy set a new folder only for these methylation sites and classified to the **Coding Region**. <br>
#### methylation site distribution
motif_methylationType.methylationLevel_result.txt

| column name | Description |
| --- | --- |
| Strains | Strain name |
| Methylation | Methylation type |
| nCDS | counts of methylated bases in CDS region |
| nPROMOTER | counts of methylated bases in promoter |
| nRR | counts of methylated bases in Regulation Region |


#### methylated gene result

motif_methylationType.methylationLevel_methylationGene.txt

| column name | Description |
| --- | --- |
| Strains | Strain name |
| Methylation | Methylation type |
| Region | the region which the methylation site located, CDS or promoter |
| RRS | Regulation region start site |
| RRE | Regulation region end site |
| Methsite | the methylation site |
| Fraction | Estimate of the fraction of molecules that carry the modification |
| distance | the distance between the start site of gene orf and methylated site |
| start | gene start site |
| end | gene end site |
| strand | the gene transcription direction, +/- |
| locus tag | the gene position ID |
| gene name | gene name annotated by prokka |
| description | gene function annotation |

motif_methylationType.methylationLevel_raw_methylationGene.txt

| column name | Description |
| --- | --- |
| Strains | Strain name |
| Methylation | Methylation type |
| Region | the region which the methylation site located, CDS or promoter |
| RRS | Regulation region start site |
| RRE | Regulation region end site |
| Methsite | the methylation site |
| distance | the distance between the start site of gene orf and methylated site |
| start | gene start site |
| end | gene end site |
| strand | the gene transcription direction, +/- |
| locus tag | the gene position ID |
| gene name | gene name annotated by prokka |
| description | gene function annotation |


#### the overlap of methylation site and TF binding
motif_methylationType.methylationLevel_TF.meme.txt

| column name | Description |
| --- | --- |
| Strains | Strain name |
| Methylation |                                                                                                                                                   Methylation type                                    |
| TF binding start |  TF binding position start site |
| TF binding end |  TF binding position end site |
| FIMO score | the TF which binding with RR or CDS region in FIMO scan score |
| Region | the region which the methylation site located, CDS or promoter |
| RRS | Regulation region start site |
| RRE | Regulation region end site |
| Methsite | the methylation site |
| Fraction | Estimate of the fraction of molecules that carry the modification |
| distance | the distance between the start site of gene orf and methylated site |
| start | gene start site |
| end | gene end site |
| strand | the gene transcription direction, +/- |
| locus tag | the gene position ID |
| gene name | gene name annotated by prokka |
| description | gene function annotation |
#### whole modified base of this motif in this genome
motif_methylationType_motif.methylationLevel<br>
e.g. GATC_m6A_motif.methylation<br>
#### Explanation
Using Calcute_PPM_console_final.py for generating meme format TFs file
---------------------
### Usage
    python ./script/Calcute_PPM_console_final.py  TF.mat
Coming updates
----------------
1. Incomplete genome assembly data<br>

Copyright
-------------------------
