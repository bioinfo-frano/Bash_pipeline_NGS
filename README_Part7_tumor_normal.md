# Part VII – Matched Tumor‑Normal Somatic Analysis (Bash Pipeline)

## Introduction

In [Part IV](README_Part4_fullbash.md) we developed a single Bash pipeline for the full **tumor‑only** somatic analysis. Although tumor-only analysis helps identify somatic variants, it relies on statistical inference and therefore provides only an approximation. In contrast, a matched tumor‑normal pair comparison is **the gold‑standard approach recommended by GATK for maximum accuracy** and, thus, is preferred whenever possible.

We will use the same reference `GRCh38` files and the same sample `SRR30536566` (tumor - colorectal cancer biopsy), but now we also include its matched normal sample `SRR30536541` from the same patient (blood).
The dataset comes from the project `PRJNA1156316` (Filipino Young‑Onset Colorectal Cancer Patients).
***The normal sample allows us to subtract germline variants and greatly reduce false positives.***

The differences between both approaches when doing somatic analysis? See **Table 1**.

### Table 1: Tumor‑Only vs. Tumor‑Normal Somatic Analysis – Overview

| Feature                  | Tumor‑Only                               | Matched Tumor‑Normal |
|--------------------------|------------------------------------------|----------------------|
| **Input**                | One BAM (tumor)                          | Two BAMs (tumor + normal from same patient) |
| **Germline handling**    | Estimated using population data (gnomAD) | Directly removed using the normal sample |
| **Normal sample**        | Not available; uses a panel of normals (PON) to filter common artifacts | Matched normal serves as a perfect baseline to subtract germline variants |
| **Core strategy**        | Statistical filtering                    | Biological comparison                    |
| **Sensitivity**          | May miss some somatic variants,<br>especially low‑allelic‑fraction | Higher sensitivity and specificity |
| **Specificity**          | Prone to false positives from sequencing artifacts and germline variants | Much higher specificity; germline variants are removed by subtraction |
| **Resources required**   | PON, germline resource (e.g., gnomAD) | Normal BAM, PON recommended (still useful for artifact filtering) |
| **Reliability**          | Moderate                                 | High (gold standard)                     |

👉 In tumor-only analysis, GATK must ***infer*** what is somatic.

👉 In tumor–normal analysis, GATK can ***directly observe*** tumor-specific variants by comparison.

---

## Where the Pipelines Are Similar

- **Pre‑processing (identical)**:
    - FastQC / MultiQC
    - Trimming (Cutadapt)
    - Alignment (BWA-MEM)
    - Sorting + MarkDuplicates + MD tags
    - BAM indexing

This part is effectively running the same preprocessing pipeline twice: once for the tumor sample and once for the matched normal.

- Post‑filtering (hard thresholds on depth, allele count, VAF) can be applied to the final somatic calls in the same way.


## Where the Pipelines Are Different

### A. Mutect2

The main differences lie in the Mutect2 command (**Table 2**).

### Table 2: Differences in Variant Calling with Mutect2 between Tumor-Only and Tumor‑Normal

| Aspect | Tumor‑Only | Matched Tumor‑Normal |
|--------|------------|----------------------|
| **Mutect2 command** | Uses `--tumor-sample` only, plus `--panel-of-normals` and `--germline-resource`. | Uses both `--tumor-sample` and `--normal-sample`. The normal BAM is provided with `-I` twice (once for each sample). The PON and germline resource are still used. |
| **Orientation bias model** | Still needed; learned from tumor BAM f1r2 tar. | Still needed; learned from tumor BAM f1r2 tar. |
| **Contamination estimation** | Estimated from tumor BAM using `GetPileupSummaries` and `CalculateContamination`. | Contamination is estimated **separately** for tumor and normal (both can be contaminated). You can run `GetPileupSummaries` on both BAMs and then `CalculateContamination` with `--matched-normal` flag. |
| **Filtering** | `FilterMutectCalls` uses contamination table and orientation model. | Same, but contamination table may include both samples. The matched normal helps in filtering germline variants. |

**Mutect2 command: tumor only**

```bash
gatk Mutect2 \
  -R reference.fasta \
  -I tumor.bam \
  --tumor-sample TUMOR_SM \
  --panel-of-normals PON \
  --germline-resource GNOMAD \
  -L intervals.bed \
  --f1r2-tar-gz f1r2.tar.gz \
  -O output.vcf.gz
```

**Mutect2 command**: tumor–normal version

```bash
gatk Mutect2 \
  -R reference.fasta \
  -I tumor.bam \
  -I normal.bam \
  --tumor-sample TUMOR_SM \
  --normal-sample NORMAL_SM \
  --panel-of-normals PON \
  --germline-resource GNOMAD \
  -L intervals.bed \
  --f1r2-tar-gz f1r2.tar.gz \
  -O output.vcf.gz
```

### B. Read group requirement

**Tumor-only**

```bash
RG_SM="DMBEL-EIDR-071"  # Library Name tumor  // Run: SRR30536566
```

**Tumor-normal**: the samples must be distinguished

```bash
# Tumor
RG_SM="DMBEL-EIDR-071"

# Normal
RG_SM="DMBEL-EIDR-096"  # Library Name normal // Run:
```
👉 The SM (sample name) must match the names provided to --tumor-sample and --normal-sample in Mutect2.


### C. Contamination estimation

**Tumor-only**

```bash
GetPileupSummaries (tumor)
→ CalculateContamination
```

**Tumor-normal**

```bash
# Tumor
gatk GetPileupSummaries -I tumor.bam ...

# Normal
gatk GetPileupSummaries -I normal.bam ...

gatk CalculateContamination \
  -I tumor.table \
  --matched-normal normal.table \
  -O contamination.table
```
👉 This allows: 
- detection of contamination in **both samples**
- more accurate correction

The contamination model becomes more accurate because allele frequencies can be compared between tumor and normal.

### D. Filtering (FilterMutectCalls)

Command is the same:

```bash
gatk FilterMutectCalls ...
```
But behavior differs:

- Tumor-only → relies on statistical filters
- Tumor–normal → uses direct evidence from the normal sample to identify germline variants and artifacts

👉 Result: fewer false positives

### E. Role of external resources

| Resource               | Tumor-only | Tumor–normal              |
| ---------------------- | ---------- | ------------------------- |
| Panel of Normals (PoN) | Essential  | Recommended               |
| gnomAD                 | Essential  | Helpful but less critical |

👉 The matched normal reduces reliance on these resources **but does not fully replace them**.

---

### Summary of Advantages or why matched tumor–normal is the gold standard

| Aspect | Matched Tumor‑Normal Benefit |
|--------|------------------------------|
| **Germline subtraction** | Removes inherited variants, leaving only somatic candidates. |
| **Better artifact filtering** | The normal helps flag sequencing errors and alignment artifacts. |
| **Contamination correction** | Both samples can be assessed for cross‑sample contamination. |
| **Sensitivity** | Better detection of low VAF somatic mutations. |
| **Higher confidence** | Final call set is smaller but more reliable. |


---

## Getting the matched tumor-normal dataset

1. Go to [SRA](https://www.ncbi.nlm.nih.gov/sra)

2. Type "SRR30536566"

3. Click on "All runs". This will send you to th "SRA RUN SELECTOR"

4...

8. 

## Folder structure

This time, the folder structure will be slightly changed. The datasets `SRR30536566` (tumor) and `SRR30536541` (matched normal) will be in the folder `~/data/PRJNA1156316_tumor_normal`. 

In addition, the subfolders `aligned`, `annotation`, `logs`, `qc`, `trimmed`, `variants` must be created for `~/SRR30536566` and `~/SRR30536541`.

In terminal, move to `~/SRR30536566`

```bash
mkdir -p aligned annotation logs qc trimmed variants
```
Continue with `~/SRR30536541`

Activate environment `DNA2`

```bash
Genomics_cancer/
├── data/
│   ├── PRJNA1156316_tumor_normal/          
│   │   └── SRR30536566                                 # TUMOR          
│   │       └── raw_fastq/
│   │           └── SRR30536566_1.fastq.gz
│   │           └── SRR30536566_2.fastq.gz           
│   │   └── SRR30536541                                 # NORMAL          
│   │       └── raw_fastq/
│   │           └── SRR30536541_1.fastq.gz
│   │           └── SRR30536541_2.fastq.gz        
│   │            
│   │
├── reference/
│   └── GRCh38/
│       ├── fasta/
│       ├── intervals/     
│       └── somatic_resources/    
│   
├── scripts/    
│   └── 09_full_somatic_PRJNA1156316_tumor_normal.nf 
│
└── logs/

```


## Conclusion

The nextflow script [09_full_somatic_SRR30536566_nextflow.nf](nextflow_scripts/09_full_somatic_SRR30536566_nextflow.nf) contains a single pipeline that runs end-to-end, performing the complete somatic analysis and producing all expected files exactly as the individual Bash scripts in [Part II – Somatic analysis](README_somatic_analysis_Part2-3.md#part-ii--somatic-analysis-bash-pipelines) and the unified Bash scripts in [Part IV](README_Part4_fullbash.md) and [Part V](README_Part5_DNA2_pipeline_update.md). The Nextflow pipeline could also output the same three expected variants, confirming reproducibility.

---

Back to the top  👉 [Part VI – Nextflow: Fully Automated Somatic DNA-NGS Pipeline](#part-vi-nextflow-pipeline-fully-automated-somatic-dna-ngs-pipeline-single-nextflow-script) 

Visit the Bash script here 👉 [Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline](README_Part4_fullbash.md)

Go and see somatic NGS analysis in `DNA2` **samtools-updated** environment in 👉 [Part V: Pipeline maintenance and Environment Validation](README_Part5_DNA2_pipeline_update.md)

Jump to the first part of this tutorial 👉 [Part I – Preparation & setup](README_setup_Part1-3.md)

Go to the main page 👉 [Bash_pipeline_NGS](README.md)

