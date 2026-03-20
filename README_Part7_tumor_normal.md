# Part VII – Matched Tumor‑Normal Somatic Analysis (Bash Pipeline)

In [Part IV](README_Part4_fullbash.md) we developed a single Bash pipeline for the full tumor‑only somatic analysis. In general, although tumor-only analysis help to identify somatic variants, this type of analysis is based on inference and, thus, provides an approximation. In contrast, a matched tumor‑normal pair comparison is **the gold‑standard approach recommended by GATK for maximum accuracy** and, thus, is preferred whenever possible.

We will use the same reference `GRCh38` files and the same sample `SRR30536566` (tumor - colorectal cancer biopsy), but now we also include its matched normal sample `SRR30536541` from the same patient (blood).
The dataset comes from the project `PRJNA1156316` (Filipino Young‑Onset Colorectal Cancer Patients).
***The normal sample allows us to subtract germline variants and greatly reduce false positives.***
In more details, what are the differences between both approaches when doing somatic analysis? See **Table 1**.

### Table 1: Tumor‑Only vs. Matched Tumor‑Normal Somatic Analysis – Overview

| Feature                  | Tumor‑Only                               | Matched Tumor‑Normal |
|--------------------------|------------------------------------------|----------------------|
| **Input**                | One BAM (tumor)                          | Two BAMs (tumor + normal from same patient) |
| **Germline handling**    | Estimated using population data (gnomAD) | Directly removed using the normal sample |
| **Normal sample**        | Not available; uses a panel of normals (PON) to filter common artifacts | Matched normal serves as a perfect baseline to subtract germline variants |
| **Core strategy**        | Statistical filtering                    | Biological comparison                    |
| **Sensitivity**          | May miss some somatic variants,<br>especially low‑allelic‑fraction | Higher sensitivity and specificity |
| **Specificity**          | Prone to false positives from sequencing artifacts and germline variants | Much higher specificity; germline variants are removed by subtraction |
| **Resources required**   | PON, germline resource (e.g., gnomAD) | Normal BAM, PON optional but still recommended |
| **Reliability**          | Moderate                                 | High (gold standard)                     |

👉 In tumor-only analysis, GATK must ***guess*** what is somatic.
👉 In tumor–normal analysis, GATK can ***observe*** the difference.

## Where the Pipelines Are Similar

- **Pre‑processing (identical)**:
    - FastQC / MultiQC
    - Trimming (Cutadapt)
    - Alignment (BWA-MEM)
    - Sorting + MarkDuplicates + MD tags
    - BAM indexing

This part is like running the script twice, once for the tumor sample and the second for the normal.

- Post‑filtering (hard thresholds on depth, allele count, VAF) can be applied to the final somatic calls in the same way.


## Differences in Variant Calling with Mutect2

### A. The main differences lie in the Mutect2 command (**Table 2**), contamination estimation, and how the matched normal is leveraged to remove germline variants.

### Table 2: Differences in Variant Calling with Mutect2 between Tumor-Only and Matched Tumor‑Normal

| Aspect | Tumor‑Only | Matched Tumor‑Normal |
|--------|------------|----------------------|
| **Mutect2 command** | Uses `--tumor-sample` only, plus `--panel-of-normals` and `--germline-resource`. | Uses both `--tumor-sample` and `--normal-sample`. The normal BAM is provided with `-I` twice (once for each sample). The PON and germline resource are still used. |
| **Orientation bias model** | Still needed; learned from tumor BAM f1r2 tar. | Still needed; learned from tumor BAM f1r2 tar. |
| **Contamination estimation** | Estimated from tumor BAM using `GetPileupSummaries` and `CalculateContamination`. | Contamination is estimated **separately** for tumor and normal (both can be contaminated). You can run `GetPileupSummaries` on both BAMs and then `CalculateContamination` with `--matched-normal` flag. |
| **Filtering** | `FilterMutectCalls` uses contamination table and orientation model. | Same, but contamination table may include both samples. The matched normal helps in filtering germline variants. |

**Mutect2 command: tumor only**

```bash
gatk Mutect2 \
  -R reference (hg38) \
  -I tumor.bam \
  --tumor-sample TUMOR (Library Name)\
  --panel-of-normals PON \
  --germline-resource GNOMAD
```

**Mutect2 command: tumor–normal version**

```bash
gatk Mutect2 \
  -I tumor.bam \
  -I normal.bam \
  --tumor-sample TUMOR \
  --normal-sample NORMAL \
  --germline-resource GNOMAD \
  --panel-of-normals PON
```

### B. Read group requirement

**Tumor-only**

```bash
RG_SM="DMBEL-EIDR-071"  # Library Name tumor  // Run: SRR30536566
```

**Tumor-normal**, the samples must be distinguished:

```bash
# Tumor
RG_SM="DMBEL-EIDR-071"

# Normal
RG_SM="DMBEL-EIDR-096"  # Library Name normal // Run:
```
👉 If both BAMs share the same sample name → Mutect2 will fail or behave incorrectly


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

### D. Filtering (FilterMutectCalls)

Command is the same:

```bash
gatk FilterMutectCalls ...
```
But behavior differs:

- Tumor-only → relies on statistical filters
- Tumor–normal → pulls direct evidence from normal

👉 Result: fewer false positives

### E. Role of external resources

| Resource               | Tumor-only | Tumor–normal              |
| ---------------------- | ---------- | ------------------------- |
| Panel of Normals (PoN) | Essential  | Recommended               |
| gnomAD                 | Essential  | Helpful but less critical |

👉 The matched normal replaces much of their function.



### Summary of Advantages or why matched tumor–normal is the gold standard

| Aspect | Matched Tumor‑Normal Benefit |
|--------|------------------------------|
| **Germline subtraction** | Removes inherited variants, leaving only somatic candidates. |
| **Better artifact filtering** | The normal helps flag sequencing errors and alignment artifacts. |
| **Contamination correction** | Both samples can be assessed for cross‑sample contamination. |
| **Sensitivity** | Better detection of low VAF somatic mutations. |
| **Higher confidence** | Final call set is smaller but more reliable. |












## 9. Outputted files `~/Genomics_cancer/data/SRR30536566_full_nf/`

```code
Genomics_cancer/
├── data/
│   ├── SRR30536566/          
│   │   └── raw_fastq/       
│   │         ├── SRR30536566_1.fastq.gz
│   │         └── SRR30536566_2.fastq.gz
│   │
│   ├── SRR30536566_full_nf/   # Used by the nextflow script
│   │   ├── qc/
│   │   │     ├── raw
│   │   │     │     ├── multiqc_data
│   │   │     │     ├── multiqc_report.html
│   │   │     │     ├── SRR30536566_1_fastqc.html
│   │   │     │     └── SRR30536566_1_fastqc.html
│   │   │     ├── md_flagstat
│   │   │     │     ├── multiqc_data
│   │   │     │     └── multiqc_report.html
│   │   │     └── trimmed
│   │   │           ├── multiqc_data
│   │   │           ├── multiqc_report.html
│   │   │           ├── SRR30536566_full_nf_R1.trimmed.fastqc.html
│   │   │           └── SRR30536566_full_nf_R2.trimmed.fastqc.html
│   │   ├── trimmed/
│   │   │     ├── SRR30536566_full_nf_R1.trimmed.fastq.gz
│   │   │     └── SRR30536566_full_nf_R2.trimmed.fastq.gz
│   │   │
│   │   ├── logs/
│   │   │     ├── cutadapt_SRR30536566_full_nf.log
│   │   │     ├── bwa_mem.log
│   │   │     ├── markduplicates.log
│   │   │     ├── SRR30536566_full_nf.flagstat.txt
│   │   │     ├── mutect2.stderr.log
│   │   │     ├── mutect2.stdout.log
│   │   │     ├── learn_read_orientation_model.log   
│   │   │     ├── get_pileup_summaries.log 
│   │   │     ├── calculate_contamination.log
│   │   │     ├── filter_mutect_calls.log
│   │   │     └── SRR30536566.postfilter.log  
│   │   ├── aligned/
│   │   │     ├── SRR30536566_full_nf.sorted.markdup.md.bam
│   │   │     ├── SRR30536566_full_nf.sorted.markdup.md.bam.bai
│   │   │     └── SRR30536566_full_nf.markdup.metrics.txt
│   │   ├── variants/
│   │   │     ├── SRR30536566_full_nf.unfiltered.vcf.gz.stats
│   │   │     ├── SRR30536566_full_nf.f1r2.tar.gz
│   │   │     ├── SRR30536566_full_nf.unfiltered.vcf.gz
│   │   │     ├── SRR30536566_full_nf.unfiltered.vcf.gz.tbi
│   │   │     ├── SRR30536566_full_nf.read-orientation-model.tar.gz
│   │   │     ├── SRR30536566_full_nf.pileups.table
│   │   │     ├── SRR30536566_full_nf.contamination.table
│   │   │     ├── SRR30536566_full_nf.filtered.vcf.gz.filteringStats.tsv
│   │   │     ├── SRR30536566_full_nf.filtered.vcf.gz
│   │   │     ├── SRR30536566_full_nf.filtered.vcf.gz.tbi
│   │   │     ├── SRR30536566_full_nf.filtered.PASS.vcf.gz
│   │   │     ├── SRR30536566_full_nf.filtered.PASS.vcf.gz.tbi
│   │   │     ├── SRR30536566_full_nf.postfiltered.vcf.gz
│   │   │     ├── SRR30536566_full_nf.postfiltered.vcf.gz.tbi
│   │   │     └── SRR30536566_full_nf.postfilter_summary.txt
│   │   └── annotation/
│   │
├── reference/
│   └── GRCh38/
│       ├── fasta/
│       ├── intervals/     
│       └── somatic_resources/    
│   
├── scripts/    
│   └── 09_full_somatic_SRR30536566_nextflow.nf 
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

