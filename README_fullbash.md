# Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline (Single Bash Script)

## Single Bash Script Implementation

In Part II, the public SRA dataset SRR30536566 (7-gene amplicon panel) was analyzed step by step using eight consecutive Bash scripts, covering:

- Quality control

- Alignment

- BAM processing

- Somatic variant calling

- Filtering and post-filtering

- Variant and clinical annotation were performed online.

In this sense:

- Part II teaches the workflow step by step.

- Part III demonstrates visualization and validation using IGV.

- Part IV now provides a fully automated, harmonized implementation of all steps in a single script.

This implementation represents a tumor-only somatic variant calling workflow.


## The full bash script

👉 [09_full_somatic_NGS_bash_script.sh](bash_scripts/09_full_somatic_NGS_bash_script.sh)

This script runs the complete NGS analysis of the SRR30536566 gene panel in a single execution, without splitting the workflow into multiple Bash files.

## What this script implements

The script integrates:

- FASTQ quality control (FastQC, MultiQC)

- Read trimming (Cutadapt)

- Alignment (BWA-MEM)

- BAM processing (sorting, duplicate marking, MD/NM tagging)

- Somatic variant calling (Mutect2)

- Orientation bias modeling

- Contamination estimation

- Variant filtering (FilterMutectCalls)

- PASS extraction

- Biological post-filtering

- Runtime reporting


## Technical Enhancements

### Error Traceability

The script identifies the exact line where an error occurs:

```bash
set -o errtrace
trap 'echo "ERROR occurred at line $LINENO"; exit 1' ERR
```

### Runtime Measurement

The script calculates total execution time:

```bash
START_TIME=$(date +%s)
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
```

## ⚠️ Prerequisites

Before running this script, you must already have:

- A correctly initialized folder structure

- Downloaded FASTQ files

- A properly indexed reference genome

- Required somatic resources (Panel of Normals, gnomAD, etc.)

For full setup instructions, see:

👉 [Part I – Preparation & setup](README_setup_Part1-3.md)


## 📁 Recommended Folder Structure

```code
Genomics_cancer/
├── data/
│   ├── SRR30536566_full/   # Used by the unified script
│   │   ├── raw_fastq/
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── logs/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   └── SRR30536566/        # Used for step-by-step workflow (Part I & II)
│       ├── raw_fastq/
│       ├── qc/
│       ├── trimmed/
│       ├── aligned/
│       ├── variants/
│       └── annotation/
│
├── reference/
│   └── GRCh38/
│       ├── fasta/
│       ├── intervals/
│       └── somatic_resources/
│
├── bash_scripts/
└── logs/
```


## 📌 Resource Dependencies

The script relies on:

- FASTQ files located in `raw_fastq/`

- Reference genome FASTA and index files in `fasta/`

- BED interval file in `intervals/`

- Panel of Normals and gnomAD resource in `somatic_resources/`














Remember that before running this script, there should be already an established folder structure and the required files, including fastq files, reference genome, gnomad among others. For all these, I recommend you to check [PART I](README_setup_Part1-3.md). Also, I recommend you create this folder structure.


Genomics_cancer/
├── data/
│   └── SRR30536566_full	     # New folder structure to run the "09_full_somatic_NGS_bash_script.sh"    
│       └── qc/          
│       └── trimmed/          
│       └── logs/          
│       └── aligned/          
│       └── variants/              
│       └── annotation/
|
│   └── SRR30536566	    # Folder structure to run the scripts step-by-step following Part I & II of this tutorial    
│       └── raw_fastq/              
│       └── qc/          
│       └── trimmed/      
│       └── aligned/          
│       └── variants/              
│       └── annotation/
├── reference/
│   └── GRCh38    
│       └── fasta/    
│       └── intervals/              
│       └── somatic_resources/                    
├── scripts/
└── logs/

The script relies on the:
- fastq files located in "raw_fastq/"    
- reference genome and fasta files located in "fasta"
- BED file from "intervals"
- Panel of Normals and gnomad from "somatic_resources" 

