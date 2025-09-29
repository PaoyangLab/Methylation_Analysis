# Bioinformatics Analysis of DNA Methylation

This repository provides a practical, end-to-end protocol for analyzing genome‑wide DNA methylation data from **bisulfite sequencing (BS‑seq)** and **enzymatic methyl sequencing (EM‑seq)**.

![DNA_methylation.png](https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/DNA_methylation.png)

---

## Table of Contents

- [Our Pipeline](#our-pipeline)
  - [Overview](#overview)
  - [Environments & Requirements](#environments--requirements)
- [Instructions & Parameters for Each Tool](#instructions--parameters-for-each-tool)
  - [QC Tools](#qc-tools)
  - [Read Processing](#read-processing)
  - [Alignment & Methylation Calling](#alignment--methylation-calling)
  - [DMR / Downstream Analysis](#dmr--downstream-analysis)
  - [Visualization](#visualization)
- [Tutorial](#tutorial)
  - [Example Dataset](#example-dataset)
  - [Step-by-Step Guidelines](#step-by-step-guidelines)
    - [1) Processing methylomes](#1-processing-methylomes)
    - [2) DMR identification](#2-dmr-identification)
    - [3) Data visualization](#3-data-visualization)
    - [4) Post-alignment analyses](#4-post-alignment-analyses)
    - [5) (Supplementary) Alternative tools](#5-supplementary-alternative-tools)
  - [Outputs (files & figures)](#outputs-files--figures)

---

## Our Pipeline

### Overview

The pipeline covers **read QC & trimming**, **alignment**, **methylation calling**, **DMR identification**, **visualization**, and **post‑alignment analyses**.

![overall_pipeline.png](https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/overall_pipeline.png)

> Main outputs include CGmap files, DMR lists (all/hyper/hypo), and various figures (PCA, heatmap, metagene plots, chromosome view).

### Environments & Requirements

Make a directory for pipeline test
```bash
mkdir ~/test_methylation_analysis
cd ~/test_methylation_analysis
```

### From GitHub
- [SRA Toolkit 3.0.5](https://github.com/ncbi/sra-tools)  
- [BS-Seeker2 v2.0.8](https://github.com/BSSeeker/BSseeker2)  
- [Trim Galore v0.4.1](https://github.com/FelixKrueger/TrimGalore)
- [HOME v1.0.0](https://github.com/ListerLab/HOME)  
- [MethylC-analyzer](https://github.com/RitataLU/MethylC-analyzer)  

> For the GitHub tools, please download the repository through `clone` command and install them by following the instructions in their manuals.
```bash
git clone "https://github.com/ncbi/sra-tools"
git clone "https://github.com/BSSeeker/BSseeker2"
git clone "https://github.com/FelixKrueger/TrimGalore"
git clone "https://github.com/ListerLab/HOME"
```

> For MethylC-analyzer, Docker image is recommended to avoid environment conflict
```bash
docker pull peiyulin/methylc:V1.0
```

###  From Other Sites 
- [Bowtie2 v2.26](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)  
- [bicycle v1.8.2](http://www.sing-group.org/bicycle)  
- [IGV Desktop v2.16.0](https://igv.org/)  
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  

> For these tools, please follow the instructions in their manuals to download and install them.

---

## Instructions & Parameters for Each Tool

> This section summarizes key usage and parameters. Full, reproducible commands are in the **Tutorial** below.

### QC Tools
**FastQC** — quality report generation  
```bash
# fastqc <seqfile(s)>
fastqc sample.fastq
```

**Trim Galore** — adapter trimming (+ auto FastQC)  
```bash
# trim_galore [--fastqc_args] "[--outdir] <output_directory>" <filename(s)>
trim_galore --fastqc_args "--outdir ./qc_trimming" sample.fastq
```

### Read Processing
**BS-Seeker2: FilterReads** — PCR duplicate removal  
```bash
# FilterReads.py -i <input> -o <output>
./BSseeker2/FilterReads.py -i sample.fastq -o sample_rmdup.fastq > sample_FilterReads.log
```

### Alignment & Methylation Calling
**BS-Seeker2: build index (bowtie2)**  
```bash
# bs_seeker2-build.py -f <reference_genome> --aligner=bowtie2 -d <output>
./BSseeker2/bs_seeker2-build.py -f genome.fa --aligner=bowtie2 -d ./BS2_bt2_Index
```

**BS-Seeker2: align**  
```bash
# bs_seeker2-align.py -i <fastq> -g <ref.fa> --aligner=bowtie2 -o <out.bam> -d <index_dir>
./BSseeker2/bs_seeker2-align.py -i sample_rmdup_trimmed.fq -g genome.fa \
  --aligner=bowtie2 -o sample_align.bam -d ./BS2_bt2_Index
```

**BS-Seeker2: call methylation**  
```bash
# bs_seeker2-call_methylation.py -i <bam> -o <out_prefix> -d <index_dir/genome.fa_bowtie2>
./BSseeker2/bs_seeker2-call_methylation.py -i sample_align.bam -o sample -d ./BS2_bt2_Index/genome.fa_bowtie2
```

**Conversion rate estimation (lambda genome)**  
- Build lambda index, align reads, call methylation as above; then:  
```bash
# conversion_rate.R <CGmap.gz>
Rscript conversion_rate.R sample_lambda.CGmap.gz
```

### DMR / Downstream Analysis
**MethylC‑analyzer (Docker)** — DMRs, Heatmap/PCA, DMG, Enrichment, Metagene, Chromosome view  
```bash
# MethylC.py [command] <sample_list> <gene.gtf> <out_dir> -a <groupA> -b <groupB>
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 \
  python /MethylC-analyzer/scripts/MethylC.py DMR samples_list.txt gene.gtf /app/ -a met1 -b wt
```

**HOME** — DMRs from CGmaps  
```bash
# HOME-pairwise -t <contexts> -i <sample_list.tsv> -o <out_dir> -mc <min_Cs> [--BSSeeker2]
HOME-pairwise -t CG -i sample_file.tsv -o ./ -mc 4 --BSSeeker2
```

**bicycle** — end-to-end methylation + differential analysis  
```bash
# bicycle create-project/reference-bisulfitation/reference-index/align/analyze-methylation
bicycle create-project -p data/myproject -r data/ref_genomes -f data/reads
bicycle reference-bisulfitation -p data/myproject
bicycle reference-index -p data/myproject -t 4
bicycle align -p data/myproject -t 4
bicycle analyze-methylation -p data/myproject -n 4 -a
```

### Visualization
**IGV Desktop** — load BigWig/TDF tracks for genome browsing.  
See Tutorial steps for tips.

---

## Tutorial

### Example Dataset
We demonstrate the pipeline using ***Arabidopsis thaliana*** dataset [GSE122394](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122394) from GEO, including wild-type (wt) strains as controls and *met1* mutant strains in which DNA methyltransferase 1 (MET1) functions primarily to maintain CG methylation. Each group (wild-type, met1 mutant) contains **three biological replicates**.  

|Name|SRA|Type|
|----|---|----|
|wt_r1|SRR8180314|Wild type|
|wt_r2|SRR8180315|Wild type|
|wt_r3|SRR8180316|Wild type|
|met1_r1|SRR8180322|MET1 mutant|
|met1_r2|SRR8180323|MET1 mutant|
|met1_r3|SRR8180324|MET1 mutant|

> To obtain the data, you can use SRA Toolkit to download the file by `prefetch` and then convert it into the FASTQ format (.fastq) for analysis by `fast-dump`.

```bash
## Download SRA data
# Usage: prefetch [options] <accessions(s)>
prefetch SRR8180314
prefetch SRR8180315
prefetch SRR8180316
prefetch SRR8180322
prefetch SRR8180323
prefetch SRR8180324

## Convert into fastq file 
#  Usage: fastq-dump [options] <accessions(s)>
#  --split-3   3-way splitting for mate-pairs. 
fastq-dump SRR8180314
fastq-dump SRR8180315
fastq-dump SRR8180316
fastq-dump SRR8180322
fastq-dump SRR8180323
fastq-dump SRR8180324
```

> Use `mv` command to rename the fastq file 

```bash
mv SRR8180314.fastq wt_r1.fastq
mv SRR8180315.fastq wt_r2.fastq
mv SRR8180316.fastq wt_r3.fastq
mv SRR8180322.fastq met1_r1.fastq
mv SRR8180323.fastq met1_r2.fastq
mv SRR8180324.fastq met1_r3.fastq
```

> Reference genome (TAIR10) can be downloaded from [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

### Step-by-Step Guidelines

### 1) Processing methylomes
Here, BS-Seeker2 is used to align reads and call methylation.

#### 1.1 Quality control and Remove duplicates
Before alignment, the methyl-seq reads should undergo quality control (QC) and trimming to remove low-quality reads and duplicate sequences generated by PCR amplification and adapter sequences. The suggested tool for QC is **FastQC**. **Trim Galore** and **BS-Seeker2** provides tools for trimming adapters and removing duplicated reads.

1.1.1 Quality control
> Fastqc generates QC report for checking read quality, the output file name will look like `wt_r1_fastqc.html` and `wt_r1_fastqc.zip`. Check the HTML file to get QC report.
```bash
# Usage: fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
#        [-c contaminant file] <seqfile(s)>
fastqc wt_r1.fastq
fastqc wt_r2.fastq
fastqc wt_r3.fastq
fastqc met1_r1.fastq
fastqc met1_r2.fastq
fastqc met1_r3.fastq
```

1.1.2 Removing duplicates
> The script `FilterReads.py` can be found from [BS-Seeker2](https://github.com/BSSeeker/BSseeker2/blob/master/FilterReads.py). The `-i` specifies the input FASTQ file like `wt_r1.fastq`, which is obtained from [Step 2](#2-download-example-data); `-o` specifies the output file `wt_r1_rmdup.fastq`; `wt_r1_FilterReads.log` records the processing log.
```bash
# Usage: FilterReads.py -i <input> -o <output>
./BSseeker2/FilterReads.py -i wt_r1.fastq -o wt_r1_rmdup.fastq > wt_r1_FilterReads.log
./BSseeker2/FilterReads.py -i wt_r2.fastq -o wt_r2_rmdup.fastq > wt_r2_FilterReads.log
./BSseeker2/FilterReads.py -i wt_r3.fastq -o wt_r3_rmdup.fastq > wt_r3_FilterReads.log
./BSseeker2/FilterReads.py -i met1_r1.fastq -o met1_r1_rmdup.fastq > met1_r1_FilterReads.log
./BSseeker2/FilterReads.py -i met1_r2.fastq -o met1_r2_rmdup.fastq > met1_r2_FilterReads.log
./BSseeker2/FilterReads.py -i met1_r3.fastq -o met1_r3_rmdup.fastq > met1_r3_FilterReads.log
```

1.1.3 Trimming adapters
> TrimGalore will automatically generate FastQC report after trimming. `--fastqc_args` set up the directory of QC report. The input file is `wt_r1_rmdup.fastq` and output files is `wt_r1_rmdup_trimmed.fq` 
```bash
## Make directory for FastQC report after trimming
mkdir qc_trimming

## Trimming
# Usage: trim_galore [--fastqc_args] "[--outdir] <output_directory>" <filename(s)>
./TrimGalore/trim_galore --fastqc_args "--outdir ./qc_trimming" wt_r1_rmdup.fastq
./TrimGalore/trim_galore --fastqc_args "--outdir ./qc_trimming" wt_r2_rmdup.fastq
./TrimGalore/trim_galore --fastqc_args "--outdir ./qc_trimming" wt_r3_rmdup.fastq
./TrimGalore/trim_galore --fastqc_args "--outdir ./qc_trimming" met1_r1_rmdup.fastq
./TrimGalore/trim_galore --fastqc_args "--outdir ./qc_trimming" met1_r2_rmdup.fastq
./TrimGalore/trim_galore --fastqc_args "--outdir ./qc_trimming" met1_r3_rmdup.fastq

## Check QC report after trimming from following directory
ls ./qc_trimming
```
#### 1.2 Alignment of methyl-seq reads
1.2.1. Use bowtie2 to create a reference genome index file (Arabidopsis thaliana TAIR10 version) for the aligner. 
> 	The `-f` specifies the FASTA file of reference genome `genome.fa`, which can be downloaded from [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html). The `-d` specifies the directory to save output index file.
```bash
## Download and untar file
wget https://s3.amazonaws.com/igenomes.illumina.com/Arabidopsis_thaliana/NCBI/TAIR10/Arabidopsis_thaliana_NCBI_TAIR10.tar.gz
tar -xzvf Arabidopsis_thaliana_NCBI_TAIR10.tar.gz

## Generate index
# Usage: bs_seeker2-build.py -f <reference_genome> --aligner=<aligner_type> -d <output>
./BSseeker2/bs_seeker2-build.py -f ./Arabidopsis_thaliana/NCBI/TAIR10/Sequence/WholeGenomeFasta/genome.fa --aligner=bowtie2 -d ./BS2_bt2_Index
```

1.2.2. Align raw reads of wild-type replicate 1 to the reference genome.
> 	The `-i` specifies input FASTQ file `wt_r1_rmdup_trimmed.fq`, which is obtained from [Step 3.1.1](#311-quality-control-and-remove-duplicates); the `-g` specifies reference genome `genome.fa` obtained from previous step; and `-o` specifies output BAM file named `wt_r1_align.bam`; `-d` specifies the index of reference genome.
```bash
# Usage: bs_seeker2-align.py -i <input_fastq> -g <reference_genome>  
# 		 --aligner=<aligner_type> -o <output_bam> -d <reference_index>
./BSseeker2/bs_seeker2-align.py -i wt_r1_rmdup_trimmed.fq -g ./Arabidopsis_thaliana/NCBI/TAIR10/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o wt_r1_align.bam -d ./BS2_bt2_Index
./BSseeker2/bs_seeker2-align.py -i wt_r2_rmdup_trimmed.fq -g ./Arabidopsis_thaliana/NCBI/TAIR10/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o wt_r2_align.bam -d ./BS2_bt2_Index
./BSseeker2/bs_seeker2-align.py -i wt_r3_rmdup_trimmed.fq -g ./Arabidopsis_thaliana/NCBI/TAIR10/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o wt_r3_align.bam -d ./BS2_bt2_Index
./BSseeker2/bs_seeker2-align.py -i met1_r1_rmdup_trimmed.fq -g ./Arabidopsis_thaliana/NCBI/TAIR10/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o met1_r1_align.bam -d ./BS2_bt2_Index
./BSseeker2/bs_seeker2-align.py -i met1_r2_rmdup_trimmed.fq -g ./Arabidopsis_thaliana/NCBI/TAIR10/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o met1_r2_align.bam -d ./BS2_bt2_Index
./BSseeker2/bs_seeker2-align.py -i met1_r3_rmdup_trimmed.fq -g ./Arabidopsis_thaliana/NCBI/TAIR10/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o met1_r3_align.bam -d ./BS2_bt2_Index
```

#### 1.3 Call methylation
1.3.1 Use methylation-calling function to calculate the methylation level. 
> 	Input BAM file `wt_r1_align.bam` is obtained from [Step 3.1.2](#312-alignment-of-methyl-seq-reads). Output file is saved as a CGmap file named `wt_r1` (the output file will be zipped into `wt_r1.CGmap.gz`). The `-d` parameter is used to specifies the index file of reference genome.
```bash
# Usage: bs_seeker2-call_methylation.py -i <input_bam> -o <output_CGmap> -d <refernce index>
./BSseeker2/bs_seeker2-call_methylation.py -i wt_r1_align.bam -o wt_r1 -d ./BS2_bt2_Index/genome.fa_bowtie2
./BSseeker2/bs_seeker2-call_methylation.py -i wt_r2_align.bam -o wt_r2 -d ./BS2_bt2_Index/genome.fa_bowtie2
./BSseeker2/bs_seeker2-call_methylation.py -i wt_r3_align.bam -o wt_r3 -d ./BS2_bt2_Index/genome.fa_bowtie2
./BSseeker2/bs_seeker2-call_methylation.py -i met1_r1_align.bam -o met1_r1 -d ./BS2_bt2_Index/genome.fa_bowtie2
./BSseeker2/bs_seeker2-call_methylation.py -i met1_r2_align.bam -o met1_r2 -d ./BS2_bt2_Index/genome.fa_bowtie2
./BSseeker2/bs_seeker2-call_methylation.py -i met1_r3_align.bam -o met1_r3 -d ./BS2_bt2_Index/genome.fa_bowtie2
```

1.3.2 View the methylation call output (CGmap). The file with each row represents a single CpG site.
```bash
zless wt_r1.CGmap.gz
```

Each CpG site contains the following information: chromosome, nucleotide on Watson strand, position, context, dinucleotide context, methylation level, number of methylated cytosines ($N_C$), and the total number of all cytosines ($N_T$ + $N_C$)

![CpGmap_table.png](https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/CpGmap_table.png)

#### 1.4 Conversion rate
In regard to methyl-seq (EM-seq and BS-seq) analysis, the estimation conversion rate, which measures how effectively bisulfite or enzyme treatment can convert unmethylated cytosines to uracil in DNA samples, is required for evaluation. By comparing the unmethylated bacteriophage lambda genome as a reference to our bisulfite/enzyme treatment genomes, the percentage of successfully converted cytosines can be estimated. It is simply calculated by dividing the number of **converted cytosines ($N_T$)** by the **total number of cytosines ($N_T$ + $N_C$)** and multiplying by 100. 

$$
\text{Conversion rate} = \frac{N_T}{N_T + N_C} \times 100\%
$$



Typically, a conversion rate of 95% or above is preferred because it shows more reliable and accurate results.

1.4.1 The first step for the conversion rate is the same as above but changes the input reference genome to the lambda phage genome.
```bash
## Download genome.fa and build index
wget https://s3.amazonaws.com/igenomes.illumina.com/Enterobacteriophage_lambda/NCBI/1993-04-28/Enterobacteriophage_lambda_NCBI_1993-04-28.tar.gz
tar -xzvf Enterobacteriophage_lambda_NCBI_1993-04-28.tar.gz

## Build index
# Usage: bs_seeker2-build.py -f <reference_genome> --aligner=<aligner_type> -d <output>
bs_seeker2-build.py -f ./Enterobacteriophage_lambda/NCBI/1993-04-28/Sequence/WholeGenomeFasta/genome.fa --aligner=bowtie2 -d ./BS2_lambda_Index
```
```bash
## Alignment
# Usage: bs_seeker2-align.py -i <input_fastq> -g <reference_genome>  
# 		 --aligner=<aligner_type> -o <output_bam> -d <reference_index>

./BSseeker2/bs_seeker2-align.py -i wt_r1_rmdup_trimmed.fq -g ./Enterobacteriophage_lambda/NCBI/1993-04-28/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o wt_r1_lambda.bam -m 3 -d ./BS2_lambda_Index
./BSseeker2/bs_seeker2-align.py -i wt_r2_rmdup_trimmed.fq -g ./Enterobacteriophage_lambda/NCBI/1993-04-28/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o wt_r2_lambda.bam -m 3 -d ./BS2_lambda_Index
./BSseeker2/bs_seeker2-align.py -i wt_r3_rmdup_trimmed.fq -g ./Enterobacteriophage_lambda/NCBI/1993-04-28/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o wt_r3_lambda.bam -m 3 -d ./BS2_lambda_Index
./BSseeker2/bs_seeker2-align.py -i met1_r1_rmdup_trimmed.fq -g ./Enterobacteriophage_lambda/NCBI/1993-04-28/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o met1_r1_lambda.bam -m 3 -d ./BS2_lambda_Index
./BSseeker2/bs_seeker2-align.py -i met1_r2_rmdup_trimmed.fq -g ./Enterobacteriophage_lambda/NCBI/1993-04-28/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o met1_r2_lambda.bam -m 3 -d ./BS2_lambda_Index
./BSseeker2/bs_seeker2-align.py -i met1_r3_rmdup_trimmed.fq -g ./Enterobacteriophage_lambda/NCBI/1993-04-28/Sequence/WholeGenomeFasta/genome.fa \
  --aligner=bowtie2 -o met1_r3_lambda.bam -m 3 -d ./BS2_lambda_Index
```
```bash
## Call methylation
# Usage: bs_seeker2-call_methylation.py -i <input_bam> -o <output_CGmap> -d <refernce index>

./BSseeker2/bs_seeker2-call_methylation.py -i wt_r1_lambda.bam -o wt_r1_lambda -d BS2_lambda_Index/genome.fa_bowtie2/
./BSseeker2/bs_seeker2-call_methylation.py -i wt_r2_lambda.bam -o wt_r2_lambda -d BS2_lambda_Index/genome.fa_bowtie2/
./BSseeker2/bs_seeker2-call_methylation.py -i wt_r3_lambda.bam -o wt_r3_lambda -d BS2_lambda_Index/genome.fa_bowtie2/
./BSseeker2/bs_seeker2-call_methylation.py -i met1_r1_lambda.bam -o met1_r1_lambda -d BS2_lambda_Index/genome.fa_bowtie2/
./BSseeker2/bs_seeker2-call_methylation.py -i met1_r2_lambda.bam -o met1_r2_lambda -d BS2_lambda_Index/genome.fa_bowtie2/
./BSseeker2/bs_seeker2-call_methylation.py -i met1_r3_lambda.bam -o met1_r3_lambda -d BS2_lambda_Index/genome.fa_bowtie2/
```

1.4.2 The conversion rate is calculated by the R script with the formula above. 
> 	The script [conversion_rate.R ](https://github.com/PaoyangLab/Methylation_Analysis/blob/main/conversion_rate.R) is included in this repository.
```bash
# Usage: conversion_rate.R <CGmap_filename>
Rscript conversion_rate.R  wt_r1_lambda.CGmap.gz
Rscript conversion_rate.R  wt_r2_lambda.CGmap.gz
Rscript conversion_rate.R  wt_r3_lambda.CGmap.gz
Rscript conversion_rate.R  met1_r1_lambda.CGmap.gz
Rscript conversion_rate.R  met1_r2_lambda.CGmap.gz
Rscript conversion_rate.R  met1_r3_lambda.CGmap.gz
```
> The output will look like (wt_r1 as example):
```bash
[ 03:13:46 AM ] Calculating bisulfite conversion rate
[ 03:13:46 AM ] Bisulfite conversion rate: 97.01493 %
```

In our example, the conversion rate for the wt_r1 methylome is 97.01%, which means that 97.01% of the unmethylated cytosines in the DNA sample have been successfully converted to uracil.

### 2) DMR identification

Here, MethylC-analyzer is selected to demonstrate how to find DMRs from the aligned methylation data output. To prevent environmental conflicts, the Docker image provided by the software is utilized
#### 2.1 Searching DMR
2.1.1 Prepare input files 
	- List CGmaps (obtained from [Step 3.1.3](#313-call-methylation)) in TXT file `samples_list.txt`, The file is tab-delimited without a header.
		> Generate `samples_list.txt`
		```bash
  		vim samples_list.txt
  		```
  		> Paste following content into `samples_list.txt`, then type `:wq` to save it.
		```tsv
		wt_r1	wt_r1.CGmap.gz	wt
		wt_r2	wt_r2.CGmap.gz	wt
		wt_r3	wt_r3.CGmap.gz	wt
  		met1_r1	met1_r1.CGmap.gz	met1
		met1_r2	met1_r2.CGmap.gz	met1
		met1_r3	met1_r3.CGmap.gz	met1
		```
	- Gene annotation GTF files `gene.gtf` can be downloaded from UCSC (https://hgdownload.soe.ucsc.edu/downloads.html), and it is a file format containing information about the genomic features of genes, such as exons, introns, coding sequences, and untranslated regions (UTRs).

2.1.2 Run MethylC-analyzer to identify the DMRs
> 	The default minimum depth for CpG sites and the number of sites within a region are both set to **4**. The default size of the DMR is **500 base pairs (bp)**. The default p-value cutoff for Student’s t-test for identifying DMRs is **P<0.05**. These arguments can be adjusted by users. The `-a` and `-b` specifies the group names.
```bash
# Usage: MethylC.py [command] <sample_list> <input_gtf_file> -a <group_a> -b <group_b>
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py DMR samples_list.txt gene.gtf /app/ -a met1 -b wt
```

The output consists of all, hyper, and hypo DMRs as text files. Here, we found **3,282 DMRs** in CG methylation between the wt and met1 groups. 

### 3) Data visualization
#### 3.1 Genome browser

![IGV_tutorial 1.png](https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/IGV_tutorial.png)

1.	Download and activate the IGV Desktop application according to the operating system. This application supports operating systems including macOS, Windows, and Linux.
2.	Select the reference genome from the dropdown list. Here, we chose A. thaliana (TAIR10) as a reference genome. Additional reference genomes can be downloaded by clicking More or can be loaded from the local path (in FASTA format).
3.	Convert the file from the WIG file to the suggested track formats, BigWig or TDF files, by running IGVtools (Click `Tools`>`Run IGVtools`).
4.	Select `File`>`Load from File` to load data into the track panel. Right-click the panel to adjust the graphic type or other settings.
5.	Use the dropdown list and search box at the top panel to select the chromosome and region shown. Click `+`/`-` on the top panel to zoom in/out. Clicking or dragging on the track of the chromosome can also adjust the region shown.
6.	Click `File`>`Save session or File`>`Save Image` to save the visualization result.

### 4) Post-alignment analyses
MethylC-analyzer provides several post-alignment analyses. Here we provide command to perform these analyses manually.

![methylC_tutorial.png](https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/methylC_tutorial.png)

> All of the following command require the same input files with [Step 3.2.1](#321-searching-dmr), including CGmap files, `sample_list.txt` and `gene.gtf`.

#### 4.1 Heatmap & PCA Analysis
```bash
# Usage: MethylC.py [command] <sample_list> <input_gtf_file> -a <group_a> -b <group_b>
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py Heatmap_PCA samples_list.txt gene.gtf /app/ -a met1 -b wt
```

The average methylation in 3 context (CG, CHG, CHH)

<img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/Average_methylation_levels.png" width="400">

PCA

<img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/PCA_CG_0.5.png" width="400">

Heatmap

<img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/Heatmap_CG_0.5.png" width="400">

#### 4.2 DMG analysis
```bash
# Usage: MethylC.py [command] <sample_list> <input_gtf_file> -a <group_a> -b <group_b>
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py DMG samples_list.txt gene.gtf /app/ -a met1 -b wt
```

Summary of dentifying Differentially Methylated Regions (DMRs) & Differentially Methylated Genes (DMGs)

<img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/Summary_DMR_DMG_numbers_CG_0.2.png" width="400">


#### 4.3 Enrichment analysis

```bash
# Usage: MethylC.py [command] <sample_list> <input_gtf_file> -a <group_a> -b <group_b>
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py Fold_Enrichment samples_list.txt gene.gtf /app/ -a met1 -b wt
```

Genomic regions fold enrichment analysis for DMRs

<img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/CG_Fold_Enrichment.png" width="400">

#### 4.4 Metagene analysis
```bash
# Usage: MethylC.py [command] <sample_list> <input_gtf_file> -a <group_a> -b <group_b>
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py Metaplot samples_list.txt gene.gtf /app/ -a met1 -b wt
```

The distribution of DNA methylation around gene body

<img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/metaplot_CG.png" width="400">

The distribution of DNA methylation difference around gene body

<img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/metaplot_delta_CG.png" width="400">

#### 4.5 Chromosome View Analysis
```bash
# Usage: MethylC.py [command] <sample_list> <input_gtf_file> -a <group_a> -b <group_b>
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py ChrView samples_list.txt gene.gtf /app/ -a met1 -b wt
```

The distribution fo DNA methylation on each chromosome

<img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/chrView_CG.png" width="800">

The distribution fo DNA methylation difference on each chromosome

<img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/chrView_delta_CG.png" width="800">


### 5) (Supplementary) Alternative tools

#### 5.1 DMR analysis by HOME
1. List CGmaps (obtained from [Step 3.1.3](#313-call-methylation)) in TXT file `sample_file.tsv`, The file is tab-delimited without a header.
> 	Generate `sample_file.tsv`
```bash
vim samples_file.tsv
```

> 	Paste following content into `samples_file.tsv`, then type `:wq` to save it.
```tsv
wt  wt1.CGmap.gz wt2.CGmap.gz wt3.CGmap.gz  
met1 met1_1.CGmap.gz met1_2.CGmap.gz met1_3.CGmap.gz
```

2. Run HOME analysis. 
> `-t` specifies methylation contexts (CG/CHG/CHH/CHN/CNN); `-i` specify input file path; `-o` specifies output directory path; `-mc` specifies minimum number of Cs in a DMR; `--BSSeeker2` indicating CGmap file from BSseeker2
```bash
# Usage: HOME-pairwise -t <contexts> -i <sample_list> -o <output_directory> -mc <min_number_of_Cs> [--BSSeeker2 source of CGmap]
HOME-pairwise -t CG -i sample_file.tsv -o ./ -mc 4 --BSSeeker2
```

#### 5.2 Methylation analysis by bicycle

1. Create a project
> `-p` specifies path to store files; `-r` specifies directory with reference genomes (put `genome.fa` from [Step3.1.2](#312-alignment-of-methyl-seq-reads) here); `-f` specifies directory with reads samples (put all `_rmdup_trimmed.fq.fq` file from [Step 3.1.1](#311-quality-control-and-remove-duplicates) here).
```bash
# Usage: bicycle [command] -p <project_path> -r <reference_genome> -f <reads_directory> 
bicycle create-project -p data/myproject -r data/ref_genomes -f data/reads
```

2. Create the Watson and Crick in-silico bisulfited reference genomes
> `-p` specifies path to store files.
```bash
# Usage: bicycle [command] -p <project_path>
bicycle reference-bisulfitation -p data/myproject
```

3. Create the bisulfited reference genome indexes
> `-p` specifies path to store files; `-t` specifies number of threads (only for bowtie2).
```bash
# Usage: bicycle [command] -p <project_path> -t <threads>
bicycle reference-index -p data/myproject -t 4
```

4. Align reads to both references
> `-p` specifies path to store files; `-t` specifies number of threads per sample and ref alignment.
```bash
# Usage: bicycle [command] -p <project_path> -t <threads_per_sample>
bicycle align -p data/myproject -t 4
```

5. Perform methylation analysis and methylcytosine calling
> `-p` specifies path to store files; `-n` specifies number of threads to analyze; `-a` ignores reads aligned to both Watson and Crick strands.
```bash
# Usage: bicycle [command] -p <project_path> -n <threads> [-a ignore double-aligned reads]
bicycle analyze-methylation -p data/myproject -n 4 -a 
```

6. Perform differential methylation analysis 
> `-p` specifies path to store files; `-t` specifies treatment-samples; `-c` specifis control-samples; `-x` specifies methylation context; `-b` specifies comma-separated (with no spaces) list of BED files to analyze at region-level.

```bash
# Usage: bicycle [command] -p <project_path> -t <treatment_sample(s)> -c <control_sample(s)> -x <context> -b <BED_file>

bicycle analyze-differential-methylation -p data/myproject_test -t met1_r1_rmdup_trimmed.fq,met1_r2_rmdup_trimmed.fq,met1_r2_rmdup_trimmed.fq -c wt_r1_rmdup_trimmed.fq,wt_r2_rmdup_trimmed.fq,wt_r2_rmdup_trimmed.fq -x CG -b TAIR10_500bp.bed
```

> The BED file for DMR analysis can generated through following bash script.
```bash
#!/bin/bash

windows=$1
Info="tair10.chrom.sizes"

bedtools makewindows -g $Info -w $windows > TAIR10_${windows}bp.bed
```
> The `tair10.chrom.sizes` is a table providing chromosom size in TAIR10 genome. The table is tab-delimited without a header:
```tsv
1       30427671
2       19698289
3       23459830
4       18585056
5       26975502
Mt      366924
Pt      154478
```

---

### Outputs (files & figures)

- **QC:** `*_fastqc.html`, logs from `FilterReads.py`, `trim_galore` outputs.  
- **Alignment:** `*_align.bam` per sample.  
- **Methylation calls:** `*.CGmap.gz` per sample; lambda `*_lambda.CGmap.gz` for conversion rate.  
- **DMRs:** text files for all/hyper/hypo (group comparison).  
- **Analyses/plots:** average methylation, PCA, heatmap, enrichment, metagene plots, chromosome view.  
- **Genome browser:** IGV session/image files.

---

## Notes
- Conversion rate (%) = `N_T / (N_T + N_C) * 100`. A value ≥95% is generally preferred.
- The example uses **A. thaliana TAIR10** as reference (download from Illumina iGenomes).
