# Part I – Preparation & setup

## Table of Contents

- [Introduction](#introduction)
- [Creating the computing environment for NGS - DNA analysis](#creating-the-computing-environment-for-ngs---dna-analysis)
    - [I. Create a specific conda environment called DNA](#i-create-a-specific-conda-environment-called-dna)
    - [II. Create the folder structure](#ii-create-the-folder-structure)
    - [III. Find & download small-sized FASTQ datasets for cancer gene panels](#iii-find--download-small-sized-fastq-datasets-for-cancer-gene-panels)
    - [IV. Download a reference human genome (GRCh38) and indexes](#iv-download-a-reference-human-genome-grch38-and-indexes)
- [Reference Genome](#reference-genome)
- [(Final) Folder structure: before starting the analysis](#final-folder-structure-before-starting-the-analysis)



## Introduction

Many students and scientists who want to learn genomic data analysis independently often rely on personal computers (e.g., laptops) with limited computational resources. These constraints can make it challenging to run secondary and tertiary NGS analyses starting from FASTQ files through variant annotation.

One alternative is to use cloud-based resources:

- **AWS EC2 instances**: Virtual machines with configurable RAM, CPU, and storage, where software dependencies can be installed via Conda. 
- **Terra.bio**: A cloud-native platform built on Google Cloud that provides pre-configured workflows, data management tools, and collaborative features specifically for biomedical research. Terra offers ready-to-use implementations of GATK and other genomics pipelines without requiring extensive infrastructure setup.
- **Cancer Genomics Cloud (CGC)**: A cloud platform powered by **Seven Bridges** that provides specialized workflows for cancer genomics, including TCGA and other cancer datasets access, pre-configured analysis pipelines, and collaborative workspaces optimized for cancer research.


However, **cloud services incur financial costs**. See below **Table 1** showing the costs for each cloud platform.

**Table 1**: Cloud platform and fees
| Platform | Website | Fee Structure & Requirements | Estimated WGS Cost (per sample) | Key Considerations & Additional Costs |
|----------|---------|-----------------------------|--------------------------------|--------------------------------------|
| **AWS** | [aws.amazon.com](https://aws.amazon.com/) | Pay-as-you-go for EC2 (compute) and S3 (storage). No platform fee. Requires AWS account and technical setup. | **$8 – $25** | • Highly instance-dependent<br>• Can drop to $3-5 with spot instances<br>• Data transfer/egress fees apply<br>• Free Tier: 12 months of limited resources |
| **Terra.bio** | [terra.bio](https://terra.bio) | Direct passthrough of **GCP (Google Cloud Platform)** rates with **no markup**. Requires active Google Cloud billing account. Optional premium features may have separate fees. | **$6 – $15** | • Costs match GCP rates exactly<br>• Preemptible VMs reduce costs by 70-80%<br>• Free Trial: $300 credit for new GCP users |
| **CGC (Cancer Genomics Cloud)** | [cancergenomicscloud.org](https://www.cancergenomicscloud.org) | Google/AWS cloud costs **plus 15-30% platform fee**. Often covered by NIH/NCI Cloud Credits for cancer researchers. Requires platform account. | **$15 – $40** | • Includes workflow management & support<br>• Optimized cancer genomics pipelines<br>• NCI Cloud Credits: Up to $10,000 annually for qualified researchers |

>**Note on Cloud Cost Estimates**: The costs above are approximate and can vary significantly based on pipeline efficiency, data volume, storage strategy, and regional pricing. Cloud platforms charge for **actual compute resource consumption per second**, storage per GB-month, and data transfer per GB, **not per sample**. These estimates are based on typical resource usage for WGS analysis and help with budgeting. All three platforms offer free tiers or credits for new users: AWS Free Tier (12 months), Google Cloud Free Trial ($300 credit), and NCI Cloud Credits for cancer researchers. For educational purposes, consider starting with small datasets or using free credits before scaling to full WGS analyses.

Therefore, for learners who prefer to work locally — typically with 8–16 GB RAM and limited disk storage (<60 GB) — it becomes essential to carefully select small sequencing datasets and design lightweight Conda environments. Equally important is the creation of a simple and efficient analysis pipeline and a well-defined folder structure.

In practice, FASTQ file sizes vary widely depending on the sequencing strategy and target region. Small targeted sequencing panels may generate FASTQ files smaller than 1 GB, whereas larger experiments such as whole-exome sequencing (WES) or whole-genome sequencing (WGS) can easily produce tens to hundreds of gigabytes per sample. For users working on local machines with limited RAM and disk space, this variability has a direct impact on dataset selection and pipeline feasibility. **Table 2** provides representative examples of FASTQ file sizes from cancer-related datasets deposited in the NCBI Sequence Read Archive (SRA).

**Table 2. Representative FASTQ file sizes from cancer-related SRA datasets**

| **Size Category** | **Estimated FASTQ Download Size (GB)** | **Sequencing Strategy** | **Example SRA Run (approx.)** | **Cancer Type / Comments** |
|------------------|----------------------------------------|-------------------------|-------------------------------|---------------------------|
| Small            | ~0.05 – 0.2                            | Targeted gene panel     | [SRX11805868](https://www.ncbi.nlm.nih.gov/sra/SRX11805868) (~0.2 GB) | Panel targeting ~95 cancer genes |
| Small            | ~0.05 – 0.1                            | Targeted gene panel     | [SRX18078826](https://www.ncbi.nlm.nih.gov/sra/SRX18078826) (~0.05 GB) | Adult breast cancer |
| Small            | ~0.06 – 0.07                           | Targeted gene panel     | [SRX18078666](https://www.ncbi.nlm.nih.gov/sra/SRX18078666) (~0.065 GB) | Adult breast cancer |
| Medium           | ~3 – 7                                 | Targeted panel / DNA-seq| [SRX28185140](https://www.ncbi.nlm.nih.gov/sra/SRX28185140) (~6 GB) | SCLC tumor sample |
| Medium / Large   | ~3 – 12                                | WES                     | [SRX29598354](https://www.ncbi.nlm.nih.gov/sra/SRX29598354) (~3.4 GB raw) | TNBC tumor |
| Large            | >50                                    | WGS                     | [SRX160842](https://www.ncbi.nlm.nih.gov/sra/SRX160842) (>>50 GB) | Acute Lymphoblastic Leukemia (ALL) Phase II |

For this reason, the examples and pipelines in this repository focus on small targeted sequencing datasets that can be processed comfortably on standard workstations or laptops.


## Creating the computing environment for NGS - DNA analysis

I.	Create conda environment

II.	Create a folder structure

III.	Find & download small-sized FASTQ datasets

IV. Download a reference human genome and indexes

V.	Bash shell scripting for NGS - DNA analysis


## Purpose of this tutorial

This repository provides a minimal, reproducible guide for running a DNA-seq (NGS) analysis pipeline locally, from FASTQ files to variant calling and annotation, using modest computational resources.

> **Note:**  
> This guide was developed and tested on macOS running on Intel processors. Users on Apple Silicon (M1/M2/…/M5) or Linux systems may need to adapt certain steps.


> **IMPORTANT – Conda prerequisites:**  
> This guide assumes that Miniconda3 is already installed on your computer. If not, please consult the official [documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/install/macos.html) or watch this [YouTube](https://www.youtube.com/watch?v=OH0E7FIHyQo) video.


When Miniconda is already installed, you should see the `(base)` environment activated in your Terminal.


## I. Create a specific conda environment called `DNA`

### a) Create the environment with these dependencies.

> **IMPORTANT**
You don't need to deactivate from `(base)` when creating a new Conda environment.

```bash
conda create -n DNA \
  -c conda-forge -c bioconda -c defaults \
  python=3.9 \
  openjdk=17 \
  perl=5.32 \
  fastqc \
  multiqc \
  cutadapt \
  bwa \
  samtools \
  picard \
  htslib \
  gatk4 \
  bcftools \
  vcftools \
  snpeff \
  ensembl-vep \
  bedtools \
  coreutils \
  pigz \
  pbzip2 \
  pandas \
  numpy \
  matplotlib \
  seaborn \
  -y
```

### b) Sanity checks - post installation

  1. List Conda environments (the new `DNA` environment should appear):

  ```bash
  conda env list
  ```

  2. Activate the new environment.

```bash
  conda activate DNA
```

  3. Verify dependencies in `DNA` and version:

```bash
  conda list
```

  gatk --version
  
  java -version
  
  bwa
  
  samtools --version
  
  vep --help
  
  snpEff -version
  

You may find the following error when running `snpEff -version`: 

```bash
Error: LinkageError
UnsupportedClassVersionError
class file version 65.0
Java runtime only recognizes up to 61.0
```
**Explanation**: 
snpEff 5.3 was compiled with Java 21
You are running Java 17 (`openjdk=17`)
Java 17 cannot run Java 21 bytecode

**Recommendation (for stability)**: 
Downgrade snpEff:

```bash
conda install -n DNA -c bioconda snpeff=5.1
```

Expected output:

```bash
SnpEff	5.1d	2022-04-19
```

## II. Create the folder structure

All FASTQ files, reference (e.g. human) genome and scripts should be located into specific folders. Below is a recommended folder structure:

```bash
Genomics_cancer/
├── reference/                 # Reference genomes and known sites
│   └── GRCh38/
│       ├── fasta/             # Reference FASTA files
│       └── known_sites/       # Known variant sites (e.g. dbSNP, Mills)
│       └── bed/               # Genomic interval files (.bed) 
│       └── somatic_resources/ # Population and somatic reference.
├── data/
│   └── SRA_ID/                # Sample-specific directory (e.g. SRR30536566)
│       ├── raw_fastq/         # Original FASTQ files
│       ├── qc/                # FastQC / MultiQC reports
│       ├── trimmed/           # Adapter- and quality-trimmed FASTQ files
│       ├── aligned/           # BAM files and indexes
│       ├── variants/          # VCF files
│       └── annotation/        # Annotated variants
├── scripts/                   # scripts
└── logs/                      # Log files from pipeline execution

bed/  
# Genomic interval files (.bed) defining target regions (e.g. gene panels,
# exons, amplicons). Used for read filtering, coverage calculation,
# and restricting variant calling to clinically relevant regions.

somatic_resources/  
# Population and somatic reference resources required for somatic variant calling,
# especially for GATK Mutect2 (e.g. Panel of Normals, gnomAD allele frequencies,
# germline resource VCFs).
```

Multiple samples can be processed by creating one directory per SRA accession under `data/`.

In Terminal, create all directories at once::

```bash
mkdir -p Genomics_cancer/{reference/GRCh38/{fasta,known_sites},data/SRR30536566/{raw_fastq,qc,trimmed,aligned,variants,annotation},scripts,logs}
```

## III. Find & download small-sized FASTQ datasets for cancer gene panels

Downloading FASTQ files directly from the SRA web interface is not recommended, because:

1. R1 and R2 reads may be merged into a single file

2. There is no guarantee that the FASTQ files represent raw reads (they may be reconstructed from alignments)

Instead, use **SRA-Tools**.

**1. Install SRA-Tools in a separate Conda environment**

```bash
conda create -n sra \
>   -c conda-forge \
>   -c bioconda \
>   python=3.10 \
>   sra-tools=3
```

Activate `(sra)` environment:

```bash
conda activate sra
```

Verify installation

```bash
fasterq-dump --version
```

Expected output:

```bash
fasterq-dump : 3.2.1			# If output is 'fasterq-dump : 2.9.6', SRA-tools won't work. Update!
```

**2. Find small FASTQ files from [SRA](https://www.ncbi.nlm.nih.gov/sra)**

> **Suggestion:**
- As an example, use these keywords: **targeted, illumina, cancer, genomic, Homo sapiens**
- As an example, you found the dataset: SRR30536566

**3. Verify dataset size and whether it contains real raw FASTQ files**.

```bash
vdb-dump --info SRR30536566
```
Expected output:

```bash
acc    : SRR30536566
path   : https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR30536566/SRR30536566
remote : https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR30536566/SRR30536566
size   : 339,709,019
type   : Database
platf  : SRA_PLATFORM_ILLUMINA
SEQ    : 3,892,036
SCHEMA : NCBI:SRA:Illumina:db#2
TIME   : 0x0000000066d7e0e7 (09/04/2024 06:24)
FMT    : sharq
FMTVER : 3.0.11
LDR    : general-loader.3.0.8
LDRVER : 3.0.8
LDRDATE: Sep 11 2023 (9/11/2023 0:0)
```
**Key fields to inspect:**
```bash
size   : 339,709,019                       # ~340 MB dataset. Compressed SRA size
SCHEMA : NCBI:SRA:Illumina:db#2            # Illumina reads forward and reverse
FMT    : sharq                             # Compressed format
```

**Table 3: Interpretation of** `vdb-dump --info`.

| **SCHEMA**                       | **FMT** | **What it really is**                              | **Suitable for FASTQ-first pipelines?** |
| -------------------------------- | ------- | -------------------------------------------------- | --------------------------------------- |
| `NCBI:SRA:GenericFastq`          | `FASTQ` | Raw FASTQ (platform-agnostic)                      | ✅ Yes                                   |
| `NCBI:SRA:GenericFastq`          | `sharq` | Raw FASTQ stored in SRA compressed format          | ✅ Yes                                   |
| `NCBI:SRA:Illumina`              | `FASTQ` | Raw FASTQ (Illumina-native representation)         | ✅ Yes                                   |
| `NCBI:SRA:Illumina`              | `sharq` | Raw Illumina FASTQ stored in compressed SRA format | ✅ Yes                                   |
| `NCBI:align:db:alignment_sorted` | `FASTQ` | FASTQ reconstructed from aligned reads             | ⚠️ Not ideal                            |
| `NCBI:align:db:alignment_sorted` | `BAM`   | Pre-aligned BAM                                    | ❌ No                                    |
| `NCBI:align:db:alignment_sorted` | `CRAM`  | Pre-aligned CRAM                                   | ❌ No                                    |

If `SCHEMA` contains `align`, the dataset is **not raw**.


**Table 4: Example of SRA dataset with raw and aligned FASTQ files.**

| Accession   | SCHEMA       | FMT   | Raw FASTQ? |
| ----------- | ------------ | ----- | ---------- |
| ERR12140864 | SRA:Illumina | sharq | ✅ Yes      | 
| SRR20701732 | align        | FASTQ | ❌ No       | 
| SRR35865210 | SRA:Illumina | sharq | ✅ Yes      | 
| SRR35529667 | SRA:Illumina | sharq | ✅ Yes      | 
| SRR32679397 | SRA:Illumina | sharq | ✅ Yes      | 

As shown in Table 3 and Table 4, both raw and aligned FASTQ datasets can be used but the `SRR20701732` is not really a raw FASTQ file.

**4. Download the SRA dataset**

- Go to working directory `~/data/SRR30536566`

```bash
fasterq-dump SRR30536566 \    # SRA Run accession ID
>   --split-files \		      	# Splits paired-end reads into two separate FASTQ files: SRR32679397_1.fastq → forward reads | SRR32679397_2.fastq → reverse reads
>   --threads 4 \		        	# Uses 4 CPU threads for parallel processing.
>   --outdir raw_fastq	  		# Specifies the output directory where the FASTQ files will be saved. Alternatively, use the full absolute path and download the dataset from any directory.
```

Expected output:

```bash
spots read      : 3,892,036
reads read      : 7,784,072
reads written   : 7,784,072
```

> **Note:**  
> The `fasterq-dump` command can be executed from any directory. However, FASTQ files should always be written to the `raw_fastq` directory using the `--outdir` option to maintain a clean and reproducible folder structure.

Then, in folder `~/raw_fastq` there should be two fastq files, each having ~1.34 GB:

SRR30536566_1.fastq

SRR30536566_2.fastq


**5. Compress the FASTQ files**

In `~/raw_fastq` directory:

```bash
gzip *.fastq
```

Expected output:

SRR30536566_1.fastq.gz

SRR30536566_2.fastq.gz


| Representation     | Approx size |
| ------------------ | ----------- |
| SRA (`sharq`)      | ~340 MB     |
| FASTQ (plain text) | ~2.6 GB (both) |
| FASTQ (gzipped)    | ~530 MB (both) |


Alternatively, compress during download with 'fastq-dump' (but it's slow, not recommended for large data)

```bash
fastq-dump SRR15506490 \
  --split-files \
  --threads 4 \
  --gzip \
  --outdir raw_fastq
```
Comparison between `fastq-dump` and `fasterq-dump`
| Tool           | Speed | Supports gzip directly | Recommendation         |
| -------------- | ----- | ---------------------- | ---------------------- |
| `fastq-dump`   | Slow  | ✅ Yes (`--gzip`)       | ❌ Avoid for large data |
| `fasterq-dump` | Fast  | ❌ No                   | ✅ Preferred tool       |



**6. Verify FASTQ integrity**

Count reads:

```bash
gzcat SRR30536566_1.fastq.gz | wc -l | awk '{print $1/4}'
```
Expected output:

3892036

```bash
gzcat SRR30536566_2.fastq.gz | wc -l | awk '{print $1/4}'
```
Expected output:

3892036

Inspect read structure:

```bash
zless SRR30536566_1.fastq.gz | head -n 8
```
Expected output:

```bash
@SRR30536566.1 K00133:507:H2N2NBBXY:7:1101:1621:1068 length=100
TNTGTTTTTCCCTCCTGTTTTTTTTTTTTTTTCCTTAAAACCTACCATTTTTTCGGCATTTTGTTTTTTTTTTTTTTTTTTTTCCTTAGATTCATAATCA
+SRR30536566.1 K00133:507:H2N2NBBXY:7:1101:1621:1068 length=100
A#AFFJJJJJJJJJJJJJJJJJJJJJJJJJJJ--<-<-77--7-7--<-7<7<77-7--<<----7A-7FF7F-<<7FAFF<--------7---7-777-
@SRR30536566.2 K00133:507:H2N2NBBXY:7:1101:2849:1086 length=100
CNCTATTCTACCGGAAGCTGGAAGCAGCCTGACAGGATGTGTGGTGCCCAAGTCTGCACAGTGAGGTGGGGAGTGAGGGCCGCAGGCAGAGGGCAAGGGA
+SRR30536566.2 K00133:507:H2N2NBBXY:7:1101:2849:1086 length=100
A#AFFFJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJFJFJJJJJJJJJJJJJJJJJJJAJJFJFJJJJJJJJJFJJJJJJJJJJJJJFJAFJJ7F<J
```

```bash
zless SRR30536566_1.fastq.gz| tail -n 8
```
Expected output:

```bash
@SRR30536566.3892035 K00133:507:H2N2NBBXY:7:2218:21339:47823 length=100
GGAGGTGGGGCCCGGTGGAGGGTGATTGGATCATGGGGGTGGATTTCTCATCAGTGGTTTAGCACTATCCCCCTTAGTGCTGTGGTCACAATAGTGAGTG
+SRR30536566.3892035 K00133:507:H2N2NBBXY:7:2218:21339:47823 length=100
AAFFFJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJFJJJ
@SRR30536566.3892036 K00133:507:H2N2NBBXY:7:2218:23267:47823 length=100
ATGGTGATTGCATCTAATGTTTTCCTGTTATAGGGCAAATAATAGTGGTGATCTGGGTAATAGTTTCTCCAAATAATGACAAGCAGAAGTATACTCTGAA
+SRR30536566.3892036 K00133:507:H2N2NBBXY:7:2218:23267:47823 length=100
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
```

```bash
zless SRR30536566_2.fastq.gz | head -n 8
```
Expected output:

```bash
@SRR30536566.1 K00133:507:H2N2NBBXY:7:1101:1621:1068 length=100
NCACAAAAACAGCAGATGAAAGAGTCCTGGAAAGCTGCCTTCAAACTACCTCAGGGTTTTCACTAACTTTTACATAGACAGCTTTGATCTTGACCAGGAA
+SRR30536566.1 K00133:507:H2N2NBBXY:7:1101:1621:1068 length=100
#-AAFJ-F77J<<J<J-<FJJFJA--7-<FJJFJF-FA---7FJJF-FJA--F7A77AFJAJJ-AJJAA-7AAJ<J-FAF<7---7F-7----7--7)FF
@SRR30536566.2 K00133:507:H2N2NBBXY:7:1101:2849:1086 length=100
NGCCTCTCCGCTGCCCTCTGCCCTGTGCCCCCTGCCCCCTGTCCCCTGTCCCTTGCCCTCTGCCTGCGGCCCTCACTCCCCACCTCACTGTGCAGACTTG
+SRR30536566.2 K00133:507:H2N2NBBXY:7:1101:2849:1086 length=100
#AAFFJFJJAJJAJJJJJJJJJJ-F<FJJJJJJ-AJJJJFJJJJJJF7FFJJJ7-FJJ<JFFJJ-FJFJJJJAJJJ7FJJJFFJAJFF7--<JF<AJAFF
```

```bash
zless SRR30536566_2.fastq.gz| tail -n 8
```
Expected output:

```bash
@SRR30536566.3892035 K00133:507:H2N2NBBXY:7:2218:21339:47823 length=100
TGGTTCTTTGGCCTATATAGACTTCTGCTTTTGGTGGGGTGGGGGGCCAGGAAGCTTCCAATCATGGCAGAAGGCAAAAAGGGAGGAGTCATCTTACACA
+SRR30536566.3892035 K00133:507:H2N2NBBXY:7:2218:21339:47823 length=100
AAAFFJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJFJJJJJJJJJJJJFJJJJJJJJJJJJJJJ<7JJAAFFFJJJJJJFF77AFJF<JJJJJ
@SRR30536566.3892036 K00133:507:H2N2NBBXY:7:2218:23267:47823 length=100
TATACTTGCCCTGATATTCTAAAACACAGAGTTTTAGTTGTTCAGAGGATAGCAACATACTTCGAGTTTTTTTCCTGATTGCTTCAGCAATTACTTGTTC
+SRR30536566.3892036 K00133:507:H2N2NBBXY:7:2218:23267:47823 length=100
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJ
```

## IV. Download a reference human genome (GRCh38) and indexes

The reference human genome is actually a bundle of files and it can be found on the GATK [website](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle), specifically in the link provided in the Resource Bundle hosted on a Google Cloud [Buckets - genomics-public-data ](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/)

In Google Cloud: ***Buckets/genomics-public-data/resources/broad/hg38/v0*** is possible to find part of the reference human genome bundle (called **hg38** (informal name) or **GRCh38** (Genome Reference Consortium human build 38)). There you can find the following files:

- Homo_sapiens_assembly38.fasta 

- Homo_sapiens_assembly38.dict 

- Homo_sapiens_assembly38.fasta.fai

The other part of the hg38 bundle can be found in Google Cloud [Buckets - gcp-public-data--broad-references](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0)

- Homo_sapiens_assembly38.fasta.64.alt

- Homo_sapiens_assembly38.fasta.64.amb

- Homo_sapiens_assembly38.fasta.64.ann

- Homo_sapiens_assembly38.fasta.64.bwt

- Homo_sapiens_assembly38.fasta.64.pac 

- Homo_sapiens_assembly38.fasta.64.sa

These are BWA index components, but now in 64-bit addressing mode.

| File      | Purpose                                   |
| --------- | ----------------------------------------- |
| `.64.amb` | Ambiguous bases (Ns)                      |
| `.64.ann` | Sequence names & lengths                  |
| `.64.bwt` | Burrows–Wheeler Transform                 |
| `.64.pac` | Packed reference sequence                 |
| `.64.sa`  | Suffix array (64-bit)                     |
| `.64.alt` | ALT contig ↔ primary contig relationships |


In Terminal, go to `/Genomics_cancer/reference/GRCh38/fasta`

**Download reference genome and indexes**: Run 👉 [0_wget_Hsapiens_assem38.sh](bash_scripts/0_wget_Hsapiens_assem38.sh)


# Reference Genome

**Genome:** GRCh38 / hg38  
**Source:** Broad Institute – GATK Resource Bundle  
**Bucket:** gcp-public-data--broad-references/hg38/v0

> **Why verify reference genome integrity?**
The reference genome FASTA and its associated index files (.fai and .dict) are foundational to all downstream NGS analyses. Any corruption or mismatch between these files can lead to alignment failures, GATK errors, or subtle coordinate inconsistencies that compromise variant calling results. Verifying file integrity using MD5 checksums ensures that the reference genome was downloaded correctly and that all associated index files correspond exactly to the same FASTA sequence, thereby guaranteeing reproducibility and analytical correctness.

MD5 checksums provide a fingerprint of a file.
If the checksum matches the expected value, the file is bit-by-bit identical to the original.

**Step 1: Go to the reference FASTA directory**

```bash
cd Genomics_cancer/reference/GRCh38/fasta
```

**Step 2: Compute MD5 checksums**

```bash
md5 Homo_sapiens_assembly38.fasta
md5 Homo_sapiens_assembly38.fasta.fai
md5 Homo_sapiens_assembly38.dict
```

Expected output:

```bash
7ff134953dcca8c8997453bbb80b6b5e  Homo_sapiens_assembly38.fasta
f76371b113734a56cde236bc0372de0a  Homo_sapiens_assembly38.fasta.fai
3884c62eb0e53fa92459ed9bff133ae6  Homo_sapiens_assembly38.dict
```


### (Final) Folder structure: before starting the analysis.

```bash
Genomics_cancer/
├── reference/                 
│   └── GRCh38/
│       ├── fasta/
│       │   ├── Homo_sapiens_assembly38.fasta
│       │   ├── Homo_sapiens_assembly38.fasta.fai
│       │   └── Homo_sapiens_assembly38.dict
│       │   └── Homo_sapiens_assembly38.fasta.64.amb     
│       │   └── Homo_sapiens_assembly38.fasta.64.ann     
│       │   └── Homo_sapiens_assembly38.fasta.64.bwt     
│       │   └── Homo_sapiens_assembly38.fasta.64.pac    
│       │   └── Homo_sapiens_assembly38.fasta.64.sa     
│       │   └── Homo_sapiens_assembly38.fasta.64.alt                  
│       └── known_sites/       
│       └── bed/               
│       └── somatic_resources/ 
├── data/
│   └── SRR30536566/                
│       ├── raw_fastq/
│       │   ├── SRR30536566_1.fastq.gz
│       │   └── SRR30536566_2.fastq.gz
│       ├── qc/                
│       ├── trimmed/           
│       ├── aligned/           
│       ├── variants/          
│       └── annotation/        
├── scripts/                   
└── logs/                      
```

If you have reached the end of **PART I**, I invite you to continue with the 👉 [Part II – Somatic analysis](README_somatic_analysis_Part2-3.md), where the DNA analysis will be explained step-by-step.

Go back to the top of 👉 [Part I – Preparation & setup](README_setup_Part1-3.md#part-i--preparation--setup)

Go to the main page 👉 [Bash_pipeline_NGS](README.md)