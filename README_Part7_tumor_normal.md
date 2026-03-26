# Part VII – Matched Tumor‑Normal Somatic Analysis (Bash Pipeline)

## Introduction

In [Part IV](README_Part4_fullbash.md) we developed a single Bash pipeline for the full **tumor‑only** somatic analysis. Although tumor-only analysis helps identify somatic variants, it relies on statistical inference and therefore provides only an approximation. In contrast, a matched tumor‑normal pair comparison is **the gold‑standard approach recommended by GATK for maximum accuracy** and, thus, is preferred whenever possible. ***The normal sample allows us to subtract germline variants and greatly reduce false positives.***

What are the differences between these two approaches in somatic analysis? See **Table 1**.

### Table 1: Tumor‑Only vs. Tumor‑Normal Somatic Analysis – Overview

| Feature                  | Tumor‑Only                               | Matched Tumor‑Normal |
|--------------------------|------------------------------------------|----------------------|
| **Input**                | One BAM (tumor)                          | Two BAMs (tumor + normal from same patient) |
| **Germline handling**    | Estimated using population data (gnomAD) | Identified and filtered using the matched normal sample |
| **Normal sample**        | Not available; uses a panel of normals (PON) to filter common artifacts | Matched normal serves as a perfect baseline to subtract germline variants |
| **Core strategy**        | Statistical filtering                    | Biological comparison                    |
| **Sensitivity**          | May miss some somatic variants,<br>especially low‑allelic‑fraction | Higher sensitivity and specificity |
| **Specificity**          | Prone to false positives from sequencing artifacts and germline variants | Much higher specificity; germline variants are removed by subtraction |
| **Resources required**   | PON, germline resource (e.g., gnomAD) | Normal BAM, PON recommended (still useful for artifact filtering) |
| **Reliability**          | Moderate                                 | High (gold standard)                     |

👉 In tumor-only analysis, GATK must ***infer*** what is somatic.

👉 In tumor–normal analysis, GATK can ***directly compare*** tumor and normal samples to identify tumor-specific variants.

---

## Where the Pipelines Are Similar

- **Pre‑processing (identical)**:
    - FastQC / MultiQC
    - Trimming (Cutadapt)
    - Alignment (BWA-MEM)
    - Sorting + MarkDuplicates + MD tags
    - BAM indexing

👉 At this stage, the tumor and normal samples are processed completely independently. They are only combined during somatic variant calling (**Mutect2**).

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


Below are simplified versions of the Mutect2 commands used in each approach:

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

👉 Mutect2 does not simply “subtract” variants found in the normal. Instead, it evaluates the likelihood that each variant is somatic given the evidence in both tumor and normal samples.

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

👉 The SM (sample name) must match the names provided to `--tumor-sample` and `--normal-sample` in Mutect2.
👉 If they do not match, Mutect2 will not correctly assign reads to tumor vs normal, leading to incorrect variant calls.

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

The contamination model becomes more accurate because information from the matched normal helps distinguish true somatic variation from contamination and germline signals.

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


### Visual workflow

A small diagram showing tumor and normal BAMs going through pre-processing independently, then merging in Mutect2

```bash
Tumor FASTQ --> Trim/Align/BAM prep \
                                   --> Mutect2 --> FilterMutectCalls --> Somatic VCF
Normal FASTQ --> Trim/Align/BAM prep /
```

---

## Matched tumor-normal: Project `PRJNA1156316`

We will use the sample `SRR30536566` (tumor - colorectal cancer biopsy) from [Part II – Somatic analysis](README_somatic_analysis_Part2-3.md), but now we also include its matched normal sample `SRR30536541` from the same patient (blood).
The dataset comes from the project `PRJNA1156316` (Filipino Young‑Onset Colorectal Cancer Patients).
- <https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1156316>


### Getting the matched tumor-normal datasets

1. Go to [SRA](https://www.ncbi.nlm.nih.gov/sra)

2. Type "SRR30536566"

>[!NOTE]
> The dataset "SRR30536566" is one of the "few" available in SRA database. Many requiere an authorised access by **dbGaP**. Details on how this dataset was found, go to chapter "**Find & download small-sized FASTQ datasets for cancer gene panels**" in  👉  [Part I - Preparation & setup](README_setup_Part1-3.md) 

3. Click on "**All runs**". This will send you to the "**SRA RUN SELECTOR**"

4. Copy the Run accession of only the matched tumor–normal pair: `SRR30536566` & `SRR30536541`

>[!NOTE]
> The datasets `SRR30536566` & `SRR30536541` must be from same patient based on "Age", "Collection_Date", "Sample Name", "sex", and "source_material_identifier"

5. Paste the Run accession, e.g. `SRR30536566` (tumor), in **[ENA](https://www.ebi.ac.uk/ena/browser/home)** "Search" pane.

6. Select the pairs "**_1**" and "**_2**" for the dataset 

7. Retrieve the `wget` by clicking in "**Get download script**"

You'll get the `wget` in a script, for example `ena-file-download-selected-files-DATE_PROJECT.sh`

8. Make the script executable with: `chmod u+x ena-file-download-selected-files-DATE_PROJECT.sh`

Repeat the download process of `wget` script for the matched normal sample.

---

### Folder structure

This time, the folder structure will be slightly different.

1. In Terminal, move to `~/Genomics_cancer/data/` and create folder `PRJNA1156316` with the subfolders `SRR30536566_tumor` and `SRR30536541_normal`:

```bash
cd ~/Genomics_cancer/data
mkdir -p PRJNA1156316/{SRR30536566_tumor,SRR30536541_normal}
```

2. In addition, the subfolder `raw_fastq` must be created for each folder `~/SRR30536566_tumor` and `~/SRR30536541_normal`.

```bash
cd ~/Genomics_cancer/data/PRJNA1156316/SRR30536566_tumor
mkdir -p raw_fastq
```
Repeat for the normal sample (`~/SRR30536541_normal`)

```bash
cd ../SRR30536541_normal
mkdir -p raw_fastq
```

3. Download the datasets: `./ena-file-download-selected-files-DATE_PROJECT.sh` in its corresponding folder:

- `~/data/PRJNA1156316/SRR30536566_tumor/raw_fastq`
- `~/data/PRJNA1156316/SRR30536541_normal/raw_fastq`


4. Activate Conda environment `DNA2`

```bash
conda activate DNA2
```

### View of folder structure for tumor-normal somatic analysis

```bash
Genomics_cancer/
├── data/
│   └── PRJNA1156316/
│       ├── SRR30536566_tumor/
│       │   └── raw_fastq/
│       │       ├── SRR30536566_1.fastq.gz
│       │       └── SRR30536566_2.fastq.gz
│       └── SRR30536541_normal
│           └── raw_fastq/
│               ├── SRR30536541_1.fastq.gz
│               └── SRR30536541_2.fastq.gz
│   
├── reference/
│   └── GRCh38/
│       ├── fasta/
│       ├── intervals/     
│       └── somatic_resources/    
│   
├── scripts/    
│   └── 09_full_somatic_DNA2_PRJNA1156316_TN.sh
│
└── logs/

```

---

### View of folder structure after running tumor-normal somatic analysis bash script 👉 [09_full_somatic_DNA2_PRJNA1156316_TN.sh](bash_scripts/09_full_somatic_DNA2_PRJNA1156316_TN.sh)

```bash
Genomics_cancer/
├── data/
│   └── PRJNA1156316/
│       ├── SRR30536566_tumor/
│       │   └── raw_fastq/
│       │   │    ├── SRR30536566_1.fastq.gz
│       │   │    └── SRR30536566_2.fastq.gz
│       │   │
│       │   ├── aligned/SRR30536566_tumor.markdup.metrics.txt, SRR30536566_tumor.sorted.markdup.md.bam, SRR30536566_tumor.sorted.markdup.md.bam.bai
│       │   ├── logs/bwa_mem.log, cutadapt_SRR30536566_tumor.log, markduplicates_SRR30536566_tumor.log, SRR30536566_tumor.flagstat.txt, SRR30536566_tumor.idxstats.txt
│       │   ├── qc/md_flagstat/, raw/, trimmed/
│       │   └── trimmed/SRR30536566_tumor_R1.trimmed.fastq.gz, SRR30536566_tumor_R2.trimmed.fastq.gz
│       │
│       └── SRR30536541_normal
│       │   └── raw_fastq/
│       │   │    ├── SRR30536541_1.fastq.gz
│       │   │    └── SRR30536541_2.fastq.gz
│       │   │
│       │   ├── aligned/SRR30536541_normal.markdup.metrics.txt, SRR30536541_normal.sorted.markdup.md.bam, SRR30536541_normal.sorted.markdup.md.bam.bai
│       │   ├── logs/bwa_mem.log, cutadapt_SRR30536541_normal.log, markduplicates_SRR30536541_normal.log, SRR30536541_normal.flagstat.txt, SRR30536541_normal.idxstats.txt
│       │   ├── qc/md_flagstat/, raw/, trimmed/
│       │   └── trimmed/SRR30536541_tumor_R1.trimmed.fastq.gz, SRR30536541_normal_R2.trimmed.fastq.gz
│       │
│       └── logs/mutect2.stderr.log, mutect2.stdout.log, PRJNA1156316_calculate_contamination.log, PRJNA1156316_filter_mutect_calls.log, PRJNA1156316_learn_read_orientation_model.log, PRJNA1156316_SRR30536566_tumor_get_pileup_summaries.log, PRJNA1156316.postfilter.log
│       │
│       └── variants/PRJNA1156316.f1r2.tar.gz, PRJNA1156316.filtered.PASS.vcf.gz, PRJNA1156316.filtered.PASS.vcf.gz.tbi, PRJNA1156316.postfiltered.vcf.gz, PRJNA1156316.postfiltered.vcf.gz.tbi, PRJNA1156316.segments.table, PRJNA1156316_normal.pileups.table, PRJNA1156316_tumor.pileups.table, PRJNA1156316.contamination.table, PRJNA1156316_tumor_normal.unfiltered.vcf.gz, PRJNA1156316.postfilter_summary.txt
│
├── reference/
│   └── GRCh38/
│       ├── fasta/
│       ├── intervals/     
│       └── somatic_resources/    
│   
└── scripts/    
    └── 09_full_somatic_DNA2_PRJNA1156316_TN.sh

```

### Results

1. Amount of post-filtered variants: **two**

```bash
cd ~/Genomics_cancer/data/PRJNA1156316/variants

cat PRJNA1156316.postfilter_summary.txt
```

Expected output:
```bash
Post-filter summary
========================
Project: PRJNA1156316
Date: Mon Mar 23 20:51:18 CET 2026

Library type: Amplicon (PCR)
Sequencing: Tumor-Normal (paired)
Post-filtering applied on tumor sample only, normal sample retained for paired analysis

Tumor sample: EIDR_55_tumor
Normal sample: EIDR_55_blood

Input VCF: /Users/Frano/Desktop/Bioinfo_2026/Genomics_cancer/data/PRJNA1156316/variants/PRJNA1156316.filtered.vcf.gz

Thresholds:
  DP >= 200
  ALT reads >= 10
  VAF >= 0.02

PASS variants before post-filtering: 2
Variants retained: 2
```

2. Type of somatic variant - genes

```bash
cd ~/Genomics_cancer/data/PRJNA1156316/variants

bcftools view -H -f PASS PRJNA1156316.postfiltered.vcf.gz
```

Expected output:
```bash
chr1	114713909	.	G	T	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=1153,1189|61,53;DP=2594;ECNT=1;ECNTH=1;GERMQ=93;MBQ=41,41;MFRL=198,157;MMQ=60,60;MPOS=25;NALOD=-3.255;NLOD=478.11;POPAF=6;ROQ=93;TLOD=318.08;AC=1;AN=4	GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/1:645,114:0.154:759:240,52:289,43:567,103:327,318,61,53	0/0:1697,0:0.0006186:1697:701,0:837,0:1590,0:826,871,0,0
chr3	179218294	.	G	A	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=929,1006|163,184;DP=2378;ECNT=1;ECNTH=1;GERMQ=93;MBQ=41,41;MFRL=191,174;MMQ=60,60;MPOS=24;NALOD=-2.988;NLOD=284.19;POPAF=6;ROQ=93;TLOD=1028.89;AC=1;AN=4	GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/1:906,347:0.28:1253:358,138:378,149:775,301:451,455,163,184	0/0:1029,0:0.001039:1029:416,0:491,0:945,0:478,551,0,0
```
The exactly same output was observed in `PRJNA1156316.filtered.vcf.gz`


**Summary of postfiltered-variants**

| Variant | Gene | Change | VAF | Depth | TLOD | Passed? | Reason |
|---------|------|--------|-----|-------|------|---------|--------|
| **chr1:114713909** | **NRAS** | G>T | 15.4% | 759 | 318.08 | ✅ PASS | Meets all thresholds |
| **chr3:179218294** | **PIK3CA** | G>A | 28.0% | 1262 | 1253 | ✅ PASS | Meets all thresholds |


---

## Conclusion

In **tumor-only** analysis (see 👉 [Part II - Variant filtering / Which variants PASS?](README_somatic_analysis_Part2-3.md)), the amount of variants after filtering were 4, after the postfiltering 3, and finally 2 after variant annotation. Similarly, the **matched tumor-normal** analysis resulted in also two variables; however, at the variant filtering level. These two variants corresponded, for both analyses: **NRAS** 181G>T / Gln61Lys and **PIK3CA** 1624G>A / Glu542Lys. Therefore, the **matched tumor-normal** analysis could detect better false positive variants and provide more reliable variant outputs.


---

👉 [Back to the top](#part-vii--matched-tumor-normal-somatic-analysis-bash-pipeline)

Go and see somatic NGS analysis in `DNA2` **samtools-updated** environment in 👉 [Part V: Pipeline maintenance and Environment Validation](README_Part5_DNA2_pipeline_update.md)

Jump to the first part of this tutorial 👉 [Part I – Preparation & setup](README_setup_Part1-3.md)

Go to the main page 👉 [Bash_pipeline_NGS](README.md)

