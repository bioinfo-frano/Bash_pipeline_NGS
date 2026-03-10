# Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline (Single Bash Script)

## Single Bash Script Implementation

In [Part II](README_somatic_analysis_Part2-3.md), the public SRA dataset `SRR30536566` (7-gene amplicon panel) was analyzed step by step using eight consecutive Bash scripts, covering:

- Quality control

- Alignment

- BAM processing

- Somatic variant calling

- Filtering and post-filtering of variants

- Variant and clinical annotation **were performed online**.

In this sense:

- Part II teaches the workflow step by step.

- Part III demonstrates visualization and validation using IGV.

- Part IV now provides a fully automated, harmonized implementation of all steps of a **tumor-only** NGS workflow analysis in a single script.


## The full bash script

👉 [09_full_somatic_NGS_bash_script.sh](bash_scripts/09_full_somatic_NGS_bash_script.sh)

This script runs the complete NGS analysis of the `SRR30536566` gene panel in a single execution, without splitting the workflow into multiple Bash files.

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


## 📁 Recommended Folder Structure for ~/SRR30536566_full

```code
Genomics_cancer/
├── data/
│   ├── SRR30536566_full/   # # NGS analysis with a single bash script `09_full_somatic_NGS_bash_script.sh` (Tutorial: Part IV). Conda `DNA` environment.
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── logs/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   └── SRR30536566/        # NGS analysis using multiple step-by-step bash scripts (Tutorial: Part I - II). Conda `DNA` environment.
│       ├── raw_fastq/
│           └── SRR30536566_1.fastq.gz
│           └── SRR30536566_2.fastq.gz
│       ├── qc/
│       ├── trimmed/
│       ├── aligned/
│       ├── variants/
│       └── annotation/
│
├── reference/
│   └── GRCh38/
│       ├── fasta/
│           └── Homo_sapiens_assembly38.fasta, .alt, .amb, .ann, .bwt, .dic, .pac, .sa, .fai
│       ├── intervals/
│           └── crc_panel_7genes_sorted.hg38.bed
│       └── somatic_resources/
│           └── 1000g_pon.hg38.vcf.gz
│           └── af-only-gnomad.hg38.vcf.gz
│
├── scripts/
└── logs/
```


## 📌 Resource Dependencies

The script relies on:

- FASTQ files located in `raw_fastq/`

- Reference genome FASTA and index files in `fasta/`

- BED interval file in `intervals/`

- Panel of Normals and gnomAD resource in `somatic_resources/`


## Running 👉 [09_full_somatic_NGS_bash_script.sh](bash_scripts/09_full_somatic_NGS_bash_script.sh): Partial view of terminal outputs

```bash
(DNA) scripts $ ls -lrth
total 80K
-rwxr--r-- 1 Frano staff  559 Jan 12 12:59  01_qc.sh
-rwxr--r-- 1 Frano staff 1.9K Jan 12 12:59  02_trim.sh
-rwxr--r-- 1 Frano staff 4.3K Jan 13 19:51 '03_align_&_bam_preprocess.sh'
-rwxr--r-- 1 Frano staff 2.3K Jan 15 08:04  04_mutect2.sh
-rwxr--r-- 1 Frano staff 1.6K Jan 15 11:41  04_make_crc_7genes_bed.sh
-rwxr--r-- 1 Frano staff 1.4K Jan 15 13:22  05_learn_read_orientation_model.sh
-rwxr--r-- 1 Frano staff  981 Jan 15 15:22  06a_get_pileup_summaries.sh
-rwxr--r-- 1 Frano staff 1.2K Jan 15 17:20  06b_calculate_contamination.sh
-rwxr--r-- 1 Frano staff 1.7K Jan 15 17:58  07_filter_mutect_calls.sh
-rwxr--r-- 1 Frano staff 3.2K Jan 25 00:27  08_postfilter.sh
-rwxr--r-- 1 Frano staff  615 Feb  3 17:52  04_for_loop_gtf.sh
-rwxr--r-- 1 Frano staff  619 Feb  3 17:56  04_for_loop_gtf_copy.sh
-rwxr--r-- 1 Frano staff  18K Feb 26 09:23  09_full_somatic_NGS_bash_script.sh
(DNA) scripts $ ./09_full_somatic_NGS_bash_script.sh
Starting FastQC...
PROJECT_ROOT: /Users/.../Genomics_cancer
SAMPLE: SRR30536566
SAMPLE_FULL: SRR30536566_full
application/gzip
application/gzip
Started analysis of SRR30536566_1.fastq.gz
Started analysis of SRR30536566_2.fastq.gz
Approx 5% complete for SRR30536566_1.fastq.gz
Approx 5% complete for SRR30536566_2.fastq.gz
...
Approx 95% complete for SRR30536566_1.fastq.gz
Approx 95% complete for SRR30536566_2.fastq.gz
Analysis complete for SRR30536566_1.fastq.gz
Analysis complete for SRR30536566_2.fastq.gz
FastQC completed successfully.

/// MultiQC v1.33

       file_search | Search path: /Users/.../Genomics_cancer/data/SRR30536566_full/qc/raw
        searching | ████████████████████████████████████████ 100% 4/4                                                   
            fastqc | Found 2 reports
     write_results | Data        : /Users/.../Genomics_cancer/data/SRR30536566_full/qc/raw/multiqc_data
     write_results | Report      : /Users/.../Genomics_cancer/data/SRR30536566_full/qc/raw/multiqc_report.html
           multiqc | MultiQC complete
MultiQC completed successfully.
Done           00:00:43     3,892,036 reads @  11.2 µs/read;   5.35 M reads/minute
Cutadapt completed successfully.
application/gzip
application/gzip
Started analysis of SRR30536566_full_R1.trimmed.fastq.gz
Started analysis of SRR30536566_full_R2.trimmed.fastq.gz
Approx 5% complete for SRR30536566_full_R1.trimmed.fastq.gz
Approx 5% complete for SRR30536566_full_R2.trimmed.fastq.gz
...
Approx 95% complete for SRR30536566_full_R1.trimmed.fastq.gz
Approx 95% complete for SRR30536566_full_R2.trimmed.fastq.gz
Analysis complete for SRR30536566_full_R1.trimmed.fastq.gz
Analysis complete for SRR30536566_full_R2.trimmed.fastq.gz
Post trimming FastQC completed successfully.

/// MultiQC v1.33

       file_search | Search path: /Users/.../Genomics_cancer/data/SRR30536566_full/qc/trimmed
        searching | ████████████████████████████████████████ 100% 4/4                                                   
            fastqc | Found 2 reports
     write_results | Data        : /Users/.../Genomics_cancer/data/SRR30536566_full/qc/trimmed/multiqc_data
     write_results | Report      : /Users/.../Genomics_cancer/data/SRR30536566_full/qc/trimmed/multiqc_report.html
           multiqc | MultiQC complete
Post trimming MultiQC completed successfully.
Checking BWA index...
BWA index found. Skipping indexing of reference genome.
Running BWA-MEM alignment...
Alignment completed.
Converting SAM to BAM...
[samopen] SAM header is present: 3366 sequences.
Sorting BAM...
Sorting completed.
Marking duplicates...
Duplicate marking completed.
Adding MD tags...
[bam_fillmd1] different NM for read 'SRR30536566.650561': 12 -> 14
[bam_fillmd1] different MD for read 'SRR30536566.650561': '18G12G0G17G5^GGGGGGA3G25' -> '18G12G0G17G5^GGGGGGA3G23N0N0'
[bam_fillmd1] different NM for read 'SRR30536566.707578': 13 -> 15
[bam_fillmd1] different MD for read 'SRR30536566.707578': '18G10G2G15G2^GGGGGGG5G1G25' -> '18G10G2G15G2^GGGGGGG5G1G23N0N0'
...
[bam_fillmd1] different NM for read 'SRR30536566.1526349': 2 -> 4
[bam_fillmd1] different MD for read 'SRR30536566.1526349': '18T1G31' -> '18T1G29N0N0'
[bam_fillmd1] different MD for read 'SRR30536566.153947': '0A36G2G3G4A14' -> '0N36G2G3G4A14'
Indexing final BAM...
BAM indexing completed.
Generating alignment statistics...
Running MultiQC...

/// MultiQC v1.33

       file_search | Search path: /Users/.../Genomics_cancer/data/SRR30536566_full/aligned
       file_search | Search path: /Users/.../Genomics_cancer/data/SRR30536566_full/logs
        searching | ████████████████████████████████████████ 100% 11/11                                                   
            picard | Found 1 MarkDuplicates reports
          samtools | Found 1 flagstat reports
          cutadapt | Found 1 reports
     write_results | Data        : /Users/.../Genomics_cancer/data/SRR30536566_full/qc/md_flagstat/multiqc_data
     write_results | Report      : /Users/.../Genomics_cancer/data/SRR30536566_full/qc/md_flagstat/multiqc_report.html
           multiqc | MultiQC complete
MultiQC of MarkDuplicates & flagstat done!
Cleaning up intermediate files...
Cleanup completed.
Alignment and BAM preprocessing completed successfully.
----------------------------------------
Alignment + BAM preprocessing completed.
Total runtime: 1319 seconds
Total runtime: 21 minutes
----------------------------------------
 
================ PREPROCESSING DONE ================

Starting somatic variant calling with Mutect2...
Running sanity checks...
All required input files found.
Running Mutect2...
Mutect2 completed successfully.
Unfiltered VCF written to:
  /Users/.../Genomics_cancer/data/SRR30536566_full/variants/SRR30536566_full.unfiltered.vcf.gz
----------------------------------------
Somatic variant calling completed.
Total runtime: 333 seconds
Total runtime: 5 minutes
----------------------------------------

Learning read orientation model starts after sanity checks...
Running sanity checks...
All required files found.
Learning read orientation model...
Tool returned:
SUCCESS
LearnReadOrientationModel completed successfully.
Orientation model written to:
  /Users/.../Genomics_cancer/data/SRR30536566_full/variants/SRR30536566_full.read-orientation-model.tar.gz
 
Starting GetPileupSummaries...
Tool returned:
SUCCESS
GetPileupSummaries completed.
 
CalculateContamination will start after sanity checks
Running sanity checks...
All required files found.
Starting CalculateContamination...
Tool returned:
SUCCESS
CalculateContamination completed.
Contamination table written to:
  /Users/.../Genomics_cancer/data/SRR30536566_full/variants/SRR30536566_full.contamination.table
 
Variant filtering will start after sanity checks
Running sanity checks...
All required files found.
Filtering Mutect2 calls...
FilterMutectCalls completed successfully.
Final filtered VCF:
  /Users/.../Genomics_cancer/data/SRR30536566_full/variants/SRR30536566_full.filtered.vcf.gz
Extracting PASS variants...
Number of PASS variants: 4
PASS-only VCF written to:
  /Users/.../Genomics_cancer/data/SRR30536566_full/variants/SRR30536566_full.filtered.PASS.vcf.gz
 
Post variant filtering will start after sanity checks
[Thu Feb 26 09:58:43 CET 2026] Starting post-filtering
Sample: SRR30536566_full
Applying post-filter thresholds:
  DP >= 200
  ALT reads (AD[1]) >= 10
  VAF >= 0.02
Indexing post-filtered VCF
Variants retained after post-filtering: 3
SRR30536566_full post-filtering completed successfully

==========================================
FULL PIPELINE COMPLETED SUCCESSFULLY
Total runtime: 1686 seconds
Total runtime: 28 minutes
Total runtime: 0 hours
==========================================
(DNA) scripts $ 
```

## Outputted files - Folder structure:

```code
Genomics_cancer/
└── data/
│    └── SRR30536566_full/
│        ├── qc/
│        │   ├── raw/
│        │       └── multiqc_report.html/multiqc_data
│        │       └── SRR30536566_[1,2]_fastqc.html/zip
│        │   ├── trimmed/
│        │       └── multiqc_report.html/multiqc_data
│        │       └── SRR30536566_full_[R1,R2].trimmed_fastqc.html/zip
│        │   └── md_flagstat/
│        │       └── multiqc_report.html/multiqc_data
│        │
│        ├── trimmed/
│        │   └── SRR30536566_full_[R1,R2].trimmed.fastq.gz
│        │
│        ├── logs/
│        │   └── cutadapt_SRR30536566_full.log
│        │   └── bwa_mem.log
│        │   └── markduplicates.log
│        │   └── SRR30536566_full.flagstat.txt
│        │   └── mutect2.stderr.log / mutect2.stdout.log
│        │   └── learn_read_orientation_model.log
│        │   └── get_pileup_summaries.log
│        │   └── calculate_contamination.log
│        │   └── filter_mutect_calls.log
│        │   └── SRR30536566_full_DNA2.postfilter.log
│        │
│        ├── aligned/
│        │   └── SRR30536566_full.sorted.markdup.md.bam
│        │   └── SRR30536566_full.sorted.markdup.md.bam.bai
│        │   └── SRR30536566_full.markdup.metrics.txt
│        │
│        ├── variants/
│        │   ├── *.contamination.table
│        │   ├── *.f1r2.tar.gz
│        │   ├── *.unfiltered.vcf.gz
│        │   ├── *.unfiltered.vcf.gz.stats
│        │   ├── *.unfiltered.vcf.gz.tbi
│        │   ├── *.read-orientation-model.tar.gz
│        │   ├── *.filtered.vcf.gz
│        │   ├── *.filtered.vcf.gz.tbi
│        │   ├── *.pileups.table
│        │   ├── *.PASS.vcf.gz
│        │   ├── *.PASS.vcf.gz.tbi
│        │   └── *.postfiltered.vcf.gz
│        │   └── *.postfiltered.vcf.gz.tbi
│        │   └── *.postfilter_summary.txt
│        │
│        └── annotation/
│
├── reference/
│   └── GRCh38/
│       ├── fasta/
│       ├── intervals/
│       └── somatic_resources/
│
├── scripts/
└── logs/
```

## Conclusion

The bash script [09_full_somatic_NGS_bash_script.sh](bash_scripts/09_full_somatic_NGS_bash_script.sh) contains a pipeline that runs smoothly and outputting all expected files as the splitted bash scripts. Most importantly, this pipeline could also output the same three expected variants.

---

Go back to the top of  👉 [Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline](README_Part4_fullbash.md.md#part-iv--bash-script-fully-automated-somatic-dna-ngs-pipeline)

Go and see somatic NGS analysis in `DNA2` **samtools-updated** environment in 👉 [Part V: Pipeline maintenance and Environment Validation](README_Part5_DNA2_pipeline_update.md)

Jump to the first part of this tutorial 👉 [Part I – Preparation & setup](README_setup_Part1-3.md)

Go to the main page 👉 [Bash_pipeline_NGS](README.md)

