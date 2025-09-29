# Bioinformatics Analysis of DNA Methylation

This repository provides a practical, end-to-end protocol for analyzing genome‑wide DNA methylation data from **bisulfite sequencing (BS‑seq)** and **enzymatic methyl sequencing (EM‑seq)**.

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

We provide two ways to set up environments: (a) native installs, or (b) Docker (recommended for MethylC‑analyzer).

**From GitHub**  
- [SRA Toolkit 3.0.5](https://github.com/ncbi/sra-tools)  
- [BS-Seeker2 v2.0.8](https://github.com/BSSeeker/BSseeker2)  
- [Trim Galore v0.4.1](https://github.com/FelixKrueger/TrimGalore)  
- [HOME v1.0.0](https://github.com/ListerLab/HOME)  
- [MethylC-analyzer](https://github.com/RitataLU/MethylC-analyzer)  

```bash
# Example: clone selected tools
git clone "https://github.com/ncbi/sra-tools"
git clone "https://github.com/BSSeeker/BSseeker2"
git clone "https://github.com/FelixKrueger/TrimGalore"
git clone "https://github.com/ListerLab/HOME"
```

**Docker image (recommended for MethylC‑analyzer)**  
```bash
docker pull peiyulin/methylc:V1.0
```

**From other sites**  
- [Bowtie2 v2.26](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)  
- [bicycle v1.8.2](http://www.sing-group.org/bicycle)  
- [IGV Desktop v2.16.0](https://igv.org/)  
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  

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

We demonstrate using ***Arabidopsis thaliana*** dataset [GSE122394](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122394) with wild‑type (wt) and *met1* mutant, each with **three biological replicates**.

| Name   | SRA       | Type       |
|--------|-----------|------------|
| wt_r1  | SRR8180314| Wild type  |
| wt_r2  | SRR8180315| Wild type  |
| wt_r3  | SRR8180316| Wild type  |
| met1_r1| SRR8180322| MET1 mutant|
| met1_r2| SRR8180323| MET1 mutant|
| met1_r3| SRR8180324| MET1 mutant|

> Download reads via `prefetch` and convert via `fastq-dump`; rename to `*_r?.fastq`.  
> Reference genome (TAIR10) can be downloaded from Illumina iGenomes.

### Step-by-Step Guidelines

#### 1) Processing methylomes
- **QC:** `fastqc <fastq>` → HTML reports in `qc_trimming` (after trimming).  
- **Remove duplicates:** `FilterReads.py -i <in.fastq> -o <out.fastq>`  
- **Trim adapters:** `trim_galore --fastqc_args "--outdir ./qc_trimming" <in.fastq>` → `*_trimmed.fq`

#### 2) DMR identification
- Prepare `samples_list.txt` (three columns: sample_id, CGmap.gz, group).  
- Run MethylC‑analyzer:
```bash
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 \
  python /MethylC-analyzer/scripts/MethylC.py DMR samples_list.txt gene.gtf /app/ -a met1 -b wt
```
- Default filters: depth/site ≥ 4; sites/region ≥ 4; region size 500bp; Student’s t‑test p<0.05 (tunable).  
- Outputs: `*_DMR_all.txt`, `*_DMR_hyper.txt`, `*_DMR_hypo.txt` (example result: **3,282 CG DMRs**).

#### 3) Data visualization
- Install **IGV Desktop**.  
- Load reference genome A. thaliana (TAIR10), convert WIG → BigWig/TDF via IGVtools.  
- `File > Load from File` to add tracks; adjust view; save session/image.

![IGV_tutorial.png](https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/IGV_tutorial.png)

#### 4) Post-alignment analyses
Run via MethylC‑analyzer:

```bash
# Heatmap & PCA
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 \
  python /MethylC-analyzer/scripts/MethylC.py Heatmap_PCA samples_list.txt gene.gtf /app/ -a met1 -b wt
# DMG
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 \
  python /MethylC-analyzer/scripts/MethylC.py DMG samples_list.txt gene.gtf /app/ -a met1 -b wt
# Enrichment
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 \
  python /MethylC-analyzer/scripts/MethylC.py Fold_Enrichment samples_list.txt gene.gtf /app/ -a met1 -b wt
# Metagene
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 \
  python /MethylC-analyzer/scripts/MethylC.py Metaplot samples_list.txt gene.gtf /app/ -a met1 -b wt
# Chromosome view
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 \
  python /MethylC-analyzer/scripts/MethylC.py ChrView samples_list.txt gene.gtf /app/ -a met1 -b wt
```

Representative figures:

- Average methylation (CG/CHG/CHH)  
  <img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/Average_methylation_levels.png" width="400">
- PCA  
  <img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/PCA_CG_0.5.png" width="400">
- Heatmap  
  <img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/Heatmap_CG_0.5.png" width="400">
- Enrichment  
  <img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/CG_Fold_Enrichment.png" width="400">
- Metagene & delta  
  <img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/metaplot_CG.png" width="400">
  <img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/metaplot_delta_CG.png" width="400">
- Chromosome view & delta  
  <img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/chrView_CG.png" width="800">
  <img src="https://github.com/PaoyangLab/Methylation_Analysis/blob/main/Figures/chrView_delta_CG.png" width="800">

#### 5) (Supplementary) Alternative tools
- **HOME** for DMRs from CGmaps.  
- **bicycle** for an alternative end‑to‑end pipeline including differential analysis.  
(See commands in **Instructions & Parameters for Each Tool**.)

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
