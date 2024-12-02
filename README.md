Bacmethy V 1.0.1
=================
Summary
------------------

Bacmethy, a customizable pipeline based on Docker, is designed to calculate the enrichment significance of methylated and un(der)methylated motifs in the regulatory and coding regions of the genome. It also identifies genes that are co-affected by DNA methylation and transcription factor binding, enabling the prediction of transcriptional regulation effects by DNA methylation. By using Bacmethy, researchers can gain a more comprehensive understanding of bacterial epigenomes.<br>
We offer three ways for utilizing Bacmethy: a web server, a Docker-based system, and a command-line tool. <br>
Please cite our work by "Liu, Ji-Hong, Yizhou Zhang, Ning Zhou, Jiale He, Jing Xu, Zhao Cai, Liang Yang, and Yang Liu. 2024. “ Bacmethy: A Novel and Convenient Tool for Investigating Bacterial DNA Methylation Pattern and Their Transcriptional Regulation Effects.” iMeta e186. https://doi.org/10.1002/imt2.186."

Installation
---------------
### 1. WEB Sever users
[Bacmethy Website](https://Bacmethy.med.sustech.edu.cn)
The website includes methylation analysis and TFs binding prediction modules of Bacmethy. <br>
No installation step is needed if using Bacmethy web server.
 
### 2. docker users
Docker image is provided for Windows or Mac users. <br>
1. Make sure you have [docker](https://docs.docker.com/engine/install/) installed. <br>
  
2. Use the command to get Bacmethy in docker:<br>
    docker pull liujihong/Bacmethy:2.0 <br>
    docker run  -t -i liujihong/Bacmethy:2.0 /bin/bash <br>
    
3. re-enter docker. <br>
    docker start <container_ID>
    docker attach <container_ID>

### 3. command-line tool users
1. install required softwares. <br>
    a. [PROKKA](https://github.com/tseemann/prokka "PROKKA Github")<br> (Note: prokka v1.12 works fine. If you want to use a higher version of prokka, you may need to manually configure the environment such as: export LD_LIBRARY_PATH = your/path/to/conda/env/lib: $LD_LIBRARY_PATH)<br>
    b. [bedtools](https://bedtools.readthedocs.io/en/latest/ "bedtools Documentation")<br>
    c. [meme](https://github.com/cinquin/MEME "meme Github")<br>
    d. r-base <br> (Note: The R packages [ggplot2](https://github.com/tidyverse/ggplot2) and [ggthemes](https://github.com/jrnold/ggthemes) need to be configured)<br>
    e. [circos](https://github.com/vigsterkr/circos) [OPTION] <br>
    f. [seqkit](https://github.com/shenwei356/seqkit) [OPTION] <br>

#### Make sure you have conda installed
    wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh  
    bash Miniconda2-latest-Linux-x86_64.sh  
#### creat a Bacmethy conda environment
    conda -n Bacmethy  
    conda activate Bacmethy  
#### conda install required softwares
    conda install -c bioconda prokka    
    conda install -c bioconda bedtools  
    conda install -c bioconda meme  
    conda install conda-forge::r-base  

2. Download the scripts from Github to your working directory. <br>
#### conda install git and clone script
    conda install git
    git clone https://github.com/LiuJih2021/Bacmethy.git
    cd ./Bacmethy/script/
      	

Inputs
-----------------
######## 
### 1. genome sequence file 
genome.fa: The referece genome sequence file (fasta file for the complete genome sequence). Note: for step 2, "Prepare motifs files," it is essential to have **complete genome** sequences obtained from the **same strain** (same reference genome). This means that the genome sequences used for this step should be derived from the exact same bacterial strains. This ensures accuracy and consistency in the analysis of motifs. 

### 2. motifs files
Nevertheless, if you need to carry out SMRTLink analysis on your own, you may download the tool and follow the instructions (https://pacbio.cn/support/software- downloads/). The tool includes the ability to assemble genomes using “Microbial Assembly Application” and detect base modifications while constructing motif‐specific models through “Base Modification Analysis Applica- tion.” When using the “Base Modification Analysis Application, it is important to add the fraction and motif requirements (‐t kineticstools_compute_methyl_fraction = true ‐t run_find_motifs = true) to obtain methylation motif and methylation fraction information accurately.<br>
1. motifs.gff: Compliant with the GFF Version3 [specification](http://www.sequenceontology.org/gff3.shtml). Each template position/strand pair whose probability value exceeds the probability value threshold appears as a row. The template position is 1-based, per the GFF specifications. The strand column refers to the strand carrying the detected modification, which is the opposite strand from those used to detect the modification. The GFF confidence column is a Phred-transformed probability value of detection.<br>
2. motifs.csv: Contains one row for each(reference,position, strand) pair that appeared in the Data Set with coverage at least x.
x defaults to 3, but is configurable with the --minCoverage option. The reference position index is 1-based for compatibility with the GFF file in the R environment. Note that this output type scales poorly and is not recommended for large genomes; the HDF5 output should perform much better in these cases.
ref: https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v10.1.pdf

### 3. PPM files [OPTION]
TF.meme: TFs (Transcription Factors) binding prediction is an optional function module of Bacmethy. To utilize this module, the TF matrix files in PPM format are needed (with a .meme suffix).

NOTE: 
Typically, you can obtain the genome.fa, motifs.gff, and motifs.csv files from a public database or through the SMRT-seq facility. The Pacbio SMRTseq facility is equipped with the pre-installed SMRTLink software. <br>
Nevertheless, if you need to carry out [SMRT-LINK](https://www.pacb.com/support/software-downloads/) analysis on your own, you may follow the instructions on the website.<br>
PacBio provide HGAP4(Hierarchical Genome Assembly Process) to generate high quality de novo assemblies of genomes, using Continuous Long Reads.<br> 
Use Base Modification Analysis Application in [SMRT-LINK](https://www.pacb.com/support/software-downloads/) or SMRT Tools to identify bacterial base modifications, and analyze the methyltransferase recognition motifs. Detection can be down using an in-silico control consisting of expected kinetic signals.<br>
When using the "Base Modification Analysis Application", it's suggested to add the fraction and motif reauirements (-t kineticstools_compute_methyl_fraction=true -t run_find_motifs=true) to obtain methylation motif and methylation fraction information accurately.<br>



Run Bacmethy
-----------------
### 1. Run Bacmethy.sh
    bash Bacmethy.sh -m motifs.gff -s motifs.csv -g genome.fasta -p PREFIX -t m6A 
#### Usage
    bash Bacmethy.sh -m <motifs.gff> -s <motifs.csv> -g <genome.fa> -p <prefix> -t <m6A, m4C or m5C> [options] -d -b -a -r -T -n -G
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
        -G      NULL    Scan the DNA methylation sites on gene Coding region (default only scan Regulation Region)
#### Notice
The genome file used for Bacmethy analysis should be the same genome reference file which used in SMRTLink methylation motif detection analysis.<br>
Additionally, there are two recommended options for better genome annotation:<br>
#### option: **-r** add standard strains reference
The **-r** option allows users to add a reference file from a standard strain in the same species. This helps improve genome annotation.
#### option: **-d** methylation level direction
The **-d** option is used for specifying the methylation level direction. There are complex and dynamic epigenetic regulations involved in DNA methylation. To ensure the authenticity of methylation events, a sequencing coverage higher than 30 is required. 
The software classifies DNA methylation events into three levels: methylation, undermethylation, and unmethylation.<br>
1. methylation: identificationQV >= 40 <br>
2. undermethylation: identificationQV >= 40 && fraction < 0.75 <br>
3. unmethylation: identifictionQV < 40 <br>
The default threshold for undermethylation in Bacmethy is set to a fraction less than 0.75 and an identificationQV greater than or equal to 40.<br>
However, users can also set user-defined thresholds using the -d parameter to change the fraction and the -i parameter to alter the identificationQV threshold.<br>

### 2. example
    Bacmethy.sh -m motifs.gff -s motif_summary.csv -g genomic.fna -p K12 -t m6A -T /Your/Path/To/TF/meme

### 3. Running Time
Bacmethy uses parallel processing to decrease running time on multicore computers.  Users can set Running CPU by parameter **-n**. 


Outputs
-----------------
structure of output files
-   methylation
    -   motif_CDS [OPTION] 
    -   motif_RR
        -   motif_methylationType.methylation (whole modified base of this motif in this genome)
        -   motif_methylationType.methylation_methylationGene.fasta (The Regulation Region in fasta format)
        -   motif_methylationType.methylation_methylationGene.txt (methylated gene result)
        -   motif_methylationType.motif.methylation_result.txt (methylation site distribution)
        -   motif_methylationType.motif.methylation_rr_locustag.txt (The corresponding locus_tag in the Regulation Region at the methylation site) 
-   undermethylation
    -   motif_CDS [OPTION] 
    -   motif_RR
-   unmethylation
    -   motif_CDS [OPTION] 
    -   motif_RR
-   PREFIX (A folder that contains all the annotated information about your sample)
-   *m6A_fraction.plot.pdf (The bar chart for fraction-frequence)
-   *m6A_plot.pdf (The scatter plot for coverage-identificationQV)
-   *m6A_motif (The raw data used for the plot)
-   *m6A_motif_basemods.txt (The clean data used for the plot)
#### Gene Features
There are 3 gene features, promoter(default: 500bp upstream region before the ATG initiation codon), CDS_RR(default: 100bp downstream region after the initiation codon), Coding Region(Whole Gene coding region). And we defined both promoter and CDS as Regulation Region(RR) in Bacmethy. <br>
##### difference between CDS_RR and Coding Region<br>
Gene transcription initiation is a complex process regulated by multiple factors. The start of the gene body can also play a role in gene transcription initiation regulation. Hence, the term **CDS_RR** is used to describe methylation sites that occur at the start of the coding sequence (CDS) region, which is involved in the regulation of gene transcription.<br>
Some users may be interested specifically in methylation sites that occur within the gene body. For this purpose, Bacmethy provides a separate folder exclusively for these methylation sites, which are classified under the category of **Coding Region**. <br>

#### methylation site distribution
motif_methylationType.motif.methylation_result.txt<br>

| column name | Description |
| --- | --- |
| Strains | Strain name |
| Methylation | Methylation type |
| nCDS | counts of methylated bases in CDS region |
| nPROMOTER | counts of methylated bases in promoter |
| nRR | counts of methylated bases in Regulation Region |


#### methylated gene result

motif_methylationType.methylation_methylationGene.txt<br>

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

#### the overlap of methylation site and TF binding

motif_methylationType.methylationLevel_TF.meme.txt<br>

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

motif_methylationType.methylation<br>
e.g. GATC_m6A_motif.methylation<br>

#### Explanation

Using Calcute_PPM_console_final.py for generating meme format TFs file<br>
---------------------
### Usage
    python ./script/Calcute_PPM_console_final.py  TF.mat

Using Circos_data_prepare.sh for generating circos configuration files and figures <br>
---------------------
### Usage (Run in the folder where the output files were generated)
    bash /PATH/TO/Circos_data_prepare.sh PREFIX /PATH/TO/SCRIPT/FOLDER/


Copyright
-------------------------
Copyright © [2024] [Ji-Hong Liu]. All rights reserved.
