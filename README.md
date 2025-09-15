# Bioinformatics Analysis of DNA Methylation

This repository provides a bioinformatics protocol for analyzing genome-wide DNA methylation data from **bisulfite sequencing (BS-seq)** and **enzymatic methyl sequencing (EM-seq)**.  

![DNA_methylation.png](https://github.com/PaoyangLab/Methyaltion_Analysis/blob/main/Figures/DNA_methylation.png)

The pipeline covers read alignment, methylation calling, DMR identification, visualization, and post-alignment analyses.  

---

## Table of Contents

- [Installation](#1-tools-installation)
- [Download Example Data](#2-download-example-data)
- [Analysis Pipeline](#3-analysis-pipeline)
	- Processing methylomes
	- DMR identification
	- Data visualization
	- Post-alignment analyses

---

## 1. Tools Installation

### From Github
- [SRA Toolkit 3.0.5](https://github.com/ncbi/sra-tools)  
- [BS-Seeker2 v2.0.8](https://github.com/BSSeeker/BSseeker2)  
- [HOME v1.0.0](https://github.com/ListerLab/HOME)  
- [MethylC-analyzer](https://github.com/RitataLU/MethylC-analyzer)  
	For MethylC-analyzer, docker image is recommended to avoid enveioment conflict
		```bash
		docker pull peiyulin/methylc:V1.0
		```
###  From Other Sites 
- [Bowtie2 v2.26](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)  
- [bicycle v1.8.2](http://www.sing-group.org/bicycle)  
- [IGV Desktop v2.16.0](https://igv.org/)  
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  

---

## 2. Download Example Data

We demonstrate the pipeline using ***Arabidopsis thaliana*** dataset (GSE122394) from GEO, including wild-type (wt) strains as controls and *met1* mutant strains in which DNA methyltransferase 1 (MET1) functions primarily to maintain CG methylation. Each group (wild-type, met1 mutant) contains **three biological replicates**.  

To obtain the data, you can use SRA Toolkit to download the file using `prefetch` and then convert it into the FASTQ format (.fastq) for analysis by `fast-dump`.
```bash
prefetch SRR8180314
fastq-dump SRR8180314
```

Reference genome (TAIR10) can be downloaded from [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

---

## 3. Analysis Pipeline

To provide useful guidance, a bioinformatics pipeline is introduced below. In the following demonstration, BS-Seeker2 is used for alignment and call methylation.

![overall_pipeline.png](https://github.com/PaoyangLab/Methyaltion_Analysis/blob/main/Figures/overall_pipeline.png)

### 3.1 Processing methylomes
#### 3.1.1 Alignments of methyl-seq read
1. Use bowtie2 to create a reference genome index file (see Note 1). for the Ar-abidopsis thaliana TAIR10 version in the aligner and save it as `BS2_bt2_Index`.
```bash
bs_seeker2-build.py -f genome.fa --aligner=bowtie2 -d ./BS2_bt2_Index
```

2. Align raw reads of wild-type replicate 1 to the reference genome using the align function and save it as a BAM file named `wt_r1_align.bam`.
```bash
bs_seeker2-align.py -i wt_r1.fastq -g genome.fa  --aligner=bowtie2 -o wt_r1_align.bam
```

#### 3.1.2 Call methylation
1. Use call methylation script to calculate the methylation level.
```bash
bs_seeker2-call_methylation.py -i wt_r1_align.bam -o wt_r1.CGmap -d ./BS2_bt2_Index/genome.fa_bowtie2
```

2. View the methylation call output (CGmap). The file with each row represents a single CpG site.
```bash
zless wt_r1.CGmap.gz
```

Each CpG site contains the following information: chromosome, nucleotide on Watson strand, position, context, dinucleotide context, methylation level, number of methylated cytosines (#C), and the total number of all cytosines (#C+T)

![CpGmap_table.png](https://github.com/PaoyangLab/Methyaltion_Analysis/blob/main/Figures/CpGmap_table.png)

#### 3.1.3 Conversion rate
1. The first step for the conversion rate is the same as above but changes the in-put reference genome to the lambda genome.
```bash
bs_seeker2-build.py -f lambda_genome.fa --aligner=bowtie2 -d ./BS2_lambda_Index

bs_seeker2-align.py -i wt_r1_rmdup.fastq -g lambda_genome.fa  --aligner=bowtie2 -o wt_r1_lambda.bam -m 3 -d BS2_lambda_Index

bs_seeker2-call_methylation.py -i wt_r1_lambda.bam -o wt_r1_lambda -d BS2_bt2_Index/genome.fa_bowtie2/
```

2. The conversion rate is calculated by the R script (see Note 3) with the formula
```bash
Rscript coversion_rate.R  wt_r1_lambda.CGmap.gz
```
```bash
# output will look like:
[ 03:13:46 AM ] Calculating bisulfite conversion rate
[ 03:13:46 AM ] Bisulfite conversion rate: 97.01493 %
```

In our example, the conversion rate for the wt_r1 methylome is 97.01%, which means that 97.01% of the unmethylated cytosines in the DNA sample have been successfully converted to uracil.

### 3.2 DMR identification
Here, MethylC-analyzer is selected to demonstrate how to find DMRs from the aligned methylation data output. To prevent environmental conflicts, the docker image provided by the software is utilized
#### 3.2.1 Searching DMR
1. Prepare input files 
	- CGmaps (samples_list.txt): obtained from Step 3.1.2
	```bash
	wt_r1.CGmap.gz
	wt_r2.CGmap.gz
	wt_r3.CGmap.gz
	met1_r1.CGmap.gz
	met1_r2.CGmap.gz
	met1_r3.CGmap.gz
	```
	- GTF files (gene.gtf): downloaded from UCSC (https://hgdownload.soe.ucsc.edu/downloads.html)
2. The default minimum depth for CpG sites and the number of sites within a region are both set to four. The default size of the DMR is 500 base pairs (bp). The default p value cutoff for Student’s t-test for identifying DMRs is 0.05. These arguments can be adjusted by users. 
```bash
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py DMR samples_list.txt gene.gtf /app/ -a met1 -b wt
```

The output consists of all, hyper, and hypo DMRs as text files. Here, we found **3,282 DMRs** in CG methylation between the wt and met1 groups. 

### 3.3	Data visualization
#### 3.3.1	Genome browser
![IGV_tutorial 1.png](https://github.com/PaoyangLab/Methyaltion_Analysis/blob/main/Figures/IGV_tutorial.png)

1.	Download and activate the IGV Desktop application according to the operat-ing system. This application supports operating systems includ-ing MacOS, Windows, and Linux.
2.	Select the reference genome from the dropdown list. Here, we chose A. thali-ana (TAIR10) as a reference genome. Additional reference genomes can be downloaded by clicking More or can be loaded from the local path (in FASTA format).
3.	Convert the file from the WIG file to the suggested track formats, BigWig or TDF files, by running IGVtools (Click Tools>Run IGVtools).
4.	Select File>Load from File to load data into the track panel. Right-click the panel to adjust the graphic type or other settings.
5.	Use the dropdown list and search box at the top panel to select the chromo-some and region shown. Click +/- on the top panel to zoom in/out. Clicking or dragging on the track of the chromosome can also adjust the region shown.
6.	Click File>Save session or File>Save Image to save the visualization result.

### 3.4	Post-alignment analyses
![methylC_tutorial.png](https://github.com/PaoyangLab/Methyaltion_Analysis/blob/main/Figures/methylC_tutorial.png)
#### 3.4.1 Enrichment analysis
Use the `Fold_Enrichment` command to generate the enrichment result. 
This module generates output files, including `CG_Fold_Enrichment.pdf` and multiple BED files, such as `CommonRegion_CG.txt.bed`. The BED format provides the information like the positions of common methylated regions across samples. The BED file can be visualized by using IGV.

```bash
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py Fold_Enrichment samples_list.txt gene.gtf /app/ -a met1 -b wt
```

DMRs exhibit a positive fold enrichment value in the IGR, suggesting a higher likelihood of DMRs being located in IGRs.

#### 3.4.2 Metagene analysis
 Use the `Metaplot` command to generate the Metaplot result. This module generates two types of metagene plots: one represents the average methylation level in two groups (metaplot_CG.pdf), and the other shows the difference 
between the two groups (metaplot_delta_CG.pdf). The former illustrates the methylation pattern along the gene body and adjacent region, while the latter directly represents the difference in distribution between wt and met1. This module also generates BigWig files (met1_r1_CG.bw) to record methylated C sites in metagene analysis, and these BigWig files can be visualized by IGV. 

```bash
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py Metaplot samples_list.txt gene.gtf /app/ -a met1 -b mt
```

In our case, the wt samples exhibit a standard CG methylation pattern [41] with a lower methylation level at the transcription start site (TSS) and transcription end site (TES). The met1 samples show a consistently low methylation level along the gene body, reflecting the dysfunction of the methyltransferase
