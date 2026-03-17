# Part V – Pipeline maintenance and Environment Validation

## Table of Contents

- [Introduction](#introduction)
- [Methodology](#methodology)
    - [1. Compare `DNA` and `DNA2` environments](#1-compare-dna-and-dna2-environments)
    - [2. Create the new conda environment `DNA2`](#2-create-the-new-conda-environment-dna2)
    - [3. Environment reproducibility (`DNA2`)](#3-environment-reproducibility-dna2)
    - [4. Verifying the integrity of the `DNA2` environment](#4-verifying-the-integrity-of-the-dna2-environment)
    - [5. Preparation of folder structure and update of Bash pipeline](#5-preparation-of-folder-structure-and-update-of-bash-pipeline)
    - [6. Run the pipeline in `DNA2`](#6-run-the-pipeline-in-dna2)
- [Results](#results)
    - [1. Filtered VCF files comparison using **BCFtools**](#filtered-vcf-files-comparison-using-bcftools)
    - [2. Post-filtered VCF files comparison using **BCFtools**](#post-filtered-vcf-files-comparison-using-bcftools)
- [Conclusion](#conclusion)

# Introduction

Maintaining or updating a bioinformatics pipeline throughout its lifecycle is an essential when working with genomic data in research or clinical contexts.

Pipeline maintenance ensures that the code written months or years ago still produces **reliable, reproducible, and scientifically valid results today**.

Because bioinformatics software evolves constantly, pipelines can gradually become outdated over time. New tool versions, updated/new reference datasets, and improved algorithms may affect analysis results. Without proper maintenance and validation, pipelines may eventually fail or produce inconsistent outputs.

Pipeline maintenance involves several types of tasks:

## 1. Tool and Dependency updates

Bioinformatics software evolves rapidly. Tools such as **GATK**, **BWA** and **samtools** are regularly updated, along with their dependencies (e.g. Python, R, or Java). For example, if a pipeline originally used **GATK 4.2** and is later updated to **GATK 4.5**, it is important to verify that the updated pipeline still produces consistent results.

Ensuring **reproducibility of the computational environment** is equally important. Package managers such as **Conda** allow the exact software environment used in the past to be recreated later. For more advanced reproducibility, container technologies such as **Docker** or **Singularity** may be used to reproduce the runtime environment.

## 2. Exceptional cases

Occasionally, a specific sequencing sample may produce an unexpected error during pipeline execution.

A robust pipeline should ideally detect these failures and allow the processing of the rest of the samples to continue uninterrupted. In practice, this may require adjusting the script to handle exceptional cases without interrupting the entire analysis run.

## 3. Reference Data Updates

Many genomic reference resources are continuously updated, including:

- dbSNP

- ClinVar

- gnomAD

- Human reference genome GRCh38.

Maintaining a pipeline may therefore require downloading updated reference datasets, indexing them, and validating compatibility with the pipeline before sharing analyses with collaborators.

## 4. Bug fixes

During routine pipeline usage, unexpected issues may arise such as:

- incorrect file paths

- software incompatibilities

- syntax changes in tool commands

- incorrect parameter usage

Fixing these issues requires updating the pipeline code while ensuring biological results remain unchanged.

## 5. Pipeline Versioning and Provenance Tracking

Tracking which pipeline version produced each result is critical. This concept is known as **pipeline versioning** and **data provenance**.

Even small changes in a pipeline may affect results. Examples include:

- Updating **samtools**

- Updating **Mutect2** in GATK

- Modifying filtering thresholds

- Updating annotation databases (gnomAD, ClinVar)

To track these changes, bioinformatics pipelines typically use **version control systems** like **Git**.

Code repositories such as **GitHub** allow developers to track every modification made to the pipeline.

In addition, pipelines often record **metadata for each run**, including:

| Metadata          | Example    |
| ----------------- | ---------- |
| Pipeline version  | v1.0       |
| Date of analysis  | 2026-03-07 |
| Reference genome  | GRCh38     |
| GATK version      | 4.6.2      |
| samtools version  | 1.22       |
| Conda environment | DNA        |


Recording this information ensures that genomic analyses remain **reproducible and traceable over time**.

## 6. Pipeline Validation

In some research institutions and diagnostic laboratories, any modification to a genomic analysis pipeline must be validated before it is used again.

Even a small change in the code or in a software dependency could potentially alter the final results. Therefore, **reference samples with known variants ("gold standard samples")** are often processed again through the updated pipeline to verify that key performance metrics such as **sensitivity** and **specificity** remain unchanged.

---

# Practical Example of Pipeline Maintenance

In this tutorial, the somatic variant analysis from a targeted NGS gene panel was performed using the Bash pipeline:

[`09_full_somatic_NGS_bash_script.sh`](bash_scripts/09_full_somatic_NGS_bash_script.sh). 

The script sequentially runs several bioinformatics tools, including the variant caller **Mutect2** from the GATK toolkit. All software runs in a Conda environment called `DNA`. 

The `DNA` environment currently contains **samtools 0.1.19**, a legacy version released before the modern **HTSlib-based samtools architecture**. Modern samtools versions (≥1.x) rely on the library **HTSlib**, , providing improved performance, new functionality, and numerous bug fixes.

An interesting question is therefore:

"***Will the number or identity of detected variants change if the pipeline is executed using a modern version of samtools?***"

---

## Methodology

To evaluate the potential impact of this update, the following approach will be used.

### 1. Compare `DNA` and `DNA2` environments

A new environment called `DNA2` will be created **with a modern version of samtools**.

**Conda environment comparison**

| Dependency | Version (`DNA`) | Dependency | Version (`DNA2`) |
| ---------- | --------------- | ---------- | ---------------- |
| Python     | 3.9             | Python     | 3.9              |
| OpenJDK    | 17.0.17         | OpenJDK    | 17.0.17          |
| Perl       | 5.32.1          | Perl       | 5.32.1           |
| samtools   | 0.1.19          | **samtools** | **1.22.1**     |
| HTSlib     | 1.22.1          | HTSlib     | 1.22.1           |
| GATK       | 4.6.2           | GATK       | 4.6.2            |
| BWA        | 0.7.19          | BWA        | 0.7.19           |
| bcftools   | 1.22            | bcftools   | 1.22             |
| Picard     | 3.4.0           | Picard     | 3.4.0            |

In this experiment, the primary tool intentionally modified is **samtools**, which is updated from version **0.1.19** to **1.22.1**.

All other core tools remain unchanged to isolate the effect of this update.

---

### 2. Create the new conda environment `DNA2`. 

The installation of `DNA2` should be performed from the `base` environment. All Conda environments are created from `base`, but remain isolated from each other.

**I. Go to `base` environment**:

```base
conda deactivate
```
You should see `(base)`

**II. Configure Conda channels (recommended)**

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

Then check:

```bash
conda config --show channels
```

Expected output:

```bash
channels:
  - conda-forge
  - bioconda
  - defaults
```

This configuration enforces a **strict channel priority**, ensuring that packages from `conda-forge` and `bioconda` are preferred over the default Anaconda channel. This reduces the risk of dependency conflicts in bioinformatics environments.

> [!CAUTION]
> During environment creation, dependency conflicts may occur when incompatible versions of tools are requested simultaneously. In such cases, removing strict version constraints often allows Conda to automatically resolve compatible versions. In genomics pipelines this frequently occurs with legacy Perl-based tools such as **Ensembl-VEP**. For this reason, complex pipelines often split tools into multiple environments to avoid dependency conflicts.

**III. Create the environment**:

```bash
conda create -n DNA2 \
  -c conda-forge -c bioconda -c defaults \
  python=3.9.23 \
  openjdk=17.0.17 \
  perl=5.32.1 \
  fastqc=0.12.1 \
  multiqc=1.33 \
  cutadapt=5.2 \
  bwa=0.7.19 \
  samtools=1.22.1 \
  picard=3.4.0 \
  htslib=1.22.1 \
  gatk4=4.6.2.0 \
  bcftools=1.22 \
  vcftools=0.1.17 \
  snpeff=5.1 \
  ensembl-vep=115.2 \
  bedtools=2.31.1 \
  coreutils=9.5 \
  pigz=2.8 \
  pbzip2=1.1.13 \
  pandas=2.3.1 \
  numpy=2.0.2 \
  matplotlib=3.9.4 \
  seaborn=0.13.2 \
  -y
```

This environment will likely produce dependency conflicts during installation. This conflict occurs because **Ensembl-VEP depends on `perl-bio-samtools`**, which in turn depends on the **legacy samtools API (<0.2)**.

Therefore, installing `ensembl-vep` forces Conda to install an outdated **samtools** version, preventing the installation of modern samtools releases (≥1.x).

```bash
ensembl-vep
  └─ perl-bioperl
       └─ perl-bio-samtools
            └─ samtools (<0.2 legacy API)
```

Since we require **samtools 1.22.1**, the dependency `ensembl-vep` must be removed and let Conda decide, which `Perl` version should be installed automatically.

The installation of `DNA2` can therefore be repeated with the following modifications.

**IV. Recreate the environment without Ensembl-VEP**

```bash
conda create -n DNA2 \
  -c conda-forge -c bioconda -c defaults \
  python=3.9.23 \
  openjdk=17.0.17 \
  perl \
  fastqc=0.12.1 \
  multiqc=1.33 \
  cutadapt=5.2 \
  bwa=0.7.19 \
  samtools=1.22.1 \
  picard=3.4.0 \
  htslib=1.22.1 \
  gatk4=4.6.2.0 \
  bcftools=1.22 \
  vcftools \
  snpeff \
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

Expected output:

```bash
Solving environment: done
...
Executing transaction: done
```

**Optional**: Instead of Conda, use **mamba**, which resolves dependencies much faster.

Check if mamba is installed:

```bash
which mamba
mamba --version
```

Then create the environment:

```bash
mamba create -n DNA2 \
  -c conda-forge -c bioconda -c defaults \
  python=3.9.23 \
  ...
```

After installation:

```bamba
conda activate DNA2
conda list
```

So now, the new `DNA2` environment will have a newer version of samtools.

| Dependency | Version (`DNA`) | Dependency   | Version (`DNA2`) |
| ---------- | --------------- | ------------ | ---------------- |
| samtools   | 0.1.19          | **samtools** | **1.22.1**       |

---

### 3. Environment reproducibility (`DNA2`)

Instead of installing environments manually, Conda environments can be exported and recreated using `.yml` or `.lock` files. Conda environments can be exported with different levels of reproducibility depending on how much information about each package is stored.

**A. Export the environment: `.yml`**

After creating `DNA2`

```bash
conda activate DNA2
conda env export --no-builds > DNA2_conda_environment.yml
```

The file `DNA2_conda_environment.yml` will store all packages and versions of `DNA2`.

For stronger reproducibility, omit `--no-builds`, which preserves exact package build identifiers.


**B. Export the environment: Create a fully reproducible `.lock` file**

If you want a perfect reproducibility, then the strongest methods are:

```bash
conda list --explicit > DNA2_conda_environment.lock
```

**C. Export the environment: for cross-platform reproducibility, use `conda-lock`** 

```bash
# First export DNA2 env without builds
conda env export --no-builds > DNA2_conda_environment.yml
# Second, create the lock DNA2 env specifically for macOS(Intel)
conda-lock lock -f DNA2_conda_environment.yml -p osx-64
# Alternatively, for macOS(Apple Silicon)
conda-lock lock -f DNA2_conda_environment.yml -p osx-arm64
# Or if you also want to increase portability to other platforms (linux-64, win-64)
conda-lock lock -f DNA2_conda_environment.yml
```

This generates:

```bash
conda-lock.yml
```

This file contains **exact package URLs and checksums** for perfect reproducibility.


**Comparison of exporting environment methods and reproducibility levels**:

| Feature                                | `conda env export`                    | `conda env export --no-builds`        | `conda list --explicit > DNA2.lock`     | `conda-lock -f DNA2_environment.yml`           |
| -------------------------------------- | ------------------------------------- | ------------------------------------- | --------------------------------------- | ---------------------------------------------- |
| **File extension**                     | `.yml`                                | `.yml`                                | `.lock` or `.txt`                       | `.lock`                                        |
| **build numbers**                      | ✔ included                            | ❌ removed                             | ✔ exact builds                          | ✔ exact builds                                 |
| **dependency solving during recreate** | minimal                               | required                              | ❌ none                                  | ❌ none (after lock created)                    |
| **reproducibility**                    | high                                  | medium                                | ⭐ perfect                               | ⭐ perfect                                      |
| **portability (other machines)**       | medium                                | high                                  | low (same platform only)                | ⭐ very high                                    |
| **cross-platform support**             | limited                               | limited                               | ❌ none                                  | ✔ macOS / Linux / clusters                     |
| **human readability**                  | ✔ easy                                | ✔ easy                                | ❌ difficult                             | medium                                         |
| **typical use**                        | sharing environment                   | sharing portable env                  | cloning exact env                       | pipelines / HPC                                |
| **used by workflow tools**             | sometimes                             | sometimes                             | rarely                                  | ✔ widely used                                  |
| **best for**                           | documentation                         | collaboration                         | exact backup                            | production workflows                           |
| **typical recreate command**           | `conda env create -f environment.yml` | `conda env create -f environment.yml` | `conda create -n DNA2 --file DNA2.lock` | `conda-lock install -n DNA2 conda-osx-64.lock` |


**D. Recreate the environment anywhere from the lock file**

To recreate the same environment in other computers and attempt to reproduce the same analysis, use these codes:

```bash
conda env create -f DNA2_conda_environment.yml
```
Expected result: a new conda environment named `DNA2` created from the YAML specification.

If the `.yml` should be directed to a different folder, for example `/envs`, then:

```bash
conda env create -f envs/DNA2_conda_environment.yml
```

If `conda-lock` was used for exporting the environment, use this code for environment recreation:

```bash
conda-lock install -n DNA2 conda-lock.yml
```
or
```bash
conda-lock -f DNA2_conda_environment.yml
```

If the environment was exported using `conda list --explicit`

```bash
conda create -n DNA2 --file env_lock.txt
```
or
```bash
conda create -n DNA2 --file env_lock.lock
```

These are **standard methods used in bioinformatics to transfer exact Conda environments**.

> [!IMPORTANT]
> By creating the environment from a **YAML** file, i.e. `conda env create -f envs/DNA2_conda_environment.yml`, Conda will **resolve dependencies again** and install compatible packages that satisfy the constraints. This means that when installing the environment through `.yml` in another system, Conda might install different builds. For example, `zlib 1.3.1` vs `zlib 1.3.1 build_1`. So the environment is reproducible **conceptually**, but not identical.
> On the other hand, by creating the environment from a **lock** file, i.e. `conda create --name DNA2 --file DNA2_lock.txt`, Conda uses an **explicit package list**. In this case, Conda **does not solve dependencies**, and it installs the **exact package builds**. So you get **bit-identical environments**. For example, `samtools-1.22.1-h96c455f_0` will be the same version, same build, same hash. The problem with this way of transferring environment is that they are **platform specific**, for example **Platform: osx-64**. This means that the installation of the environment **would fail** on:
> - Linux cluster
> - Apple Silicon
> - HPC systems
>
> YAML files are portable, shareable, and widely used in bioinformatics workflows, but they are not guaranteed to produce bit-identical environments because dependencies are solved again during installation.

---

### 4. Verifying the integrity of the `DNA2` environment

A “clean” environment means that all packages are resolved correctly with no hidden conflicts.

Before using the environment `DNA2`:

**I. Activate the new environment**:

```bash
conda activate DNA2
```

**II. Check the installation history**:

```bash
conda list --revisions
```
output:

```bash
2026-03-08 11:02:06  (rev 0)
    +_openmp_mutex-4.5 (conda-forge/osx-64)
    +_python_abi3_support-1.0 (conda-forge/noarch)
    ...
    +zstandard-0.23.0 (conda-forge/osx-64)
    +zstd-1.5.7 (conda-forge/osx-64)
```
**Meaning** of `(rev 0)`: It means that `DNA2` environment has never been modified since creation.

```bash
rev 0  -> environment created
rev 1  -> you installed a package
rev 2  -> you upgraded something
```

**III. Simulate installing critical packages to check for conflicts**:

```bash
mamba install --dry-run -n DNA2 samtools htslib bcftools
```

output:

```bash
Looking for: ['samtools', 'htslib', 'bcftools']

bioconda/noarch                                               No change
pkgs/main/osx-64                                              No change
bioconda/osx-64                                               No change
pkgs/r/osx-64                                                 No change
pkgs/r/noarch                                                 No change
pkgs/main/noarch                                              No change
conda-forge/noarch                                  25.0MB @   1.3MB/s 18.8s
conda-forge/osx-64                                  44.5MB @   2.2MB/s 20.0s

Pinned packages:
  - python 3.9.*


Transaction

  Prefix: /opt/anaconda3/envs/DNA2

  All requested packages already installed
```
**Meaning** of `All requested packages already installed`: It means that samtools, htslib, bcftools were already installed in `DNA2` environment, and there is **no dependency resolution needed**. No conflicts indicate a clean environment.


> [!Note]
> Using `--dry-run` ensures that you detect dependency conflicts **without modifying the environment**, preventing surprises later in your analysis.


---

### 5. Preparation of folder structure and Update of Bash pipeline

**I. Make a copy of the bash script `09_full_somatic_NGS_bash_script.sh`**

Go to `~/Genomics_cancer/scripts`

```bash
cp 09_full_somatic_NGS_bash_script.sh 09_full_somatic_DNA2_updated.sh
```

**II. Create a new folder called "SRR30536566_full_DNA2" with its subfolders**

Go to `~/Genomics_cancer/data`

```bash
mkdir -p SRR30536566_full_DNA2/{aligned,logs,qc,trimmed,variants,annotation}
```

Folder structure:

```code
Genomics_cancer/
├── data/
│   ├── SRR30536566_full/   # NGS analysis with a single bash script `09_full_somatic_NGS_bash_script.sh` (Tutorial: Part IV). DNA environment.
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── logs/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   └── SRR30536566/        # NGS analysis using multiple step-by-step bash scripts (Tutorial: Part I - II). DNA environment.
│   │   ├── raw_fastq/
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   └── SRR30536566_full_DNA2/        # NGS analysis using a single samtools-updated bash script (Tutorial: Part V). DNA2 environment.
│       ├── qc/
│       ├── trimmed/
│       ├── logs/
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
├── scripts/
└── logs/
```

**III Update the bash script `09_full_somatic_DNA2_updated.sh`**

Change variables in "Configuration"

```bash
SAMPLE="SRR30536566"
SAMPLE_FULL="${SAMPLE}_full"
```

To this one:

```bash
SAMPLE="SRR30536566"
SAMPLE_FULL="${SAMPLE}_full_DNA2"
```

Change samtools syntax in "Step3: Sort BAM (coordinate sort)"

From this version:

```bash
samtools sort -@ "$THREADS" \
  "$ALIGN_DIR/${SAMPLE_FULL}.bam" \
  "$ALIGN_DIR/${SAMPLE_FULL}.sorted"
```

To this one
```bash
samtools sort -@ "$THREADS" \
  -o "$ALIGN_DIR/${SAMPLE_FULL}.sorted.bam" \
  "$ALIGN_DIR/${SAMPLE_FULL}.bam"
```

**Everything else automatically updates**, this means that the pipeline in `09_full_somatic_DNA2_updated.sh` will output to `data/SRR30536566_full_DNA2/` as intended.

---

### 6. Run the pipeline in `DNA2`

**I. Activate the environment**

```bash
conda activate DNA2
```

**II. Move to the `scripts` directory and execute the pipeline: `09_full_somatic_DNA2_updated.sh`**

```bash
cd ~/Genomics_cancer/scripts
./09_full_somatic_DNA2_updated.sh
```
The pipeline will process the same sequencing dataset `SRR30536566`.

### Expected output & folder structure

```code
Genomics_cancer/
└── data/
    └── SRR30536566_full_DNA2/
        ├── qc/
        │   ├── raw/
        │       └── multiqc_report.html/multiqc_data
        │       └── SRR30536566_[1,2]_fastqc.html/zip
        │   ├── trimmed/
        │       └── multiqc_report.html/multiqc_data
        │       └── SRR30536566_full_DNA2_[R1,R2].trimmed_fastqc.html/zip
        │   └── md_flagstat/
        │       └── multiqc_report.html/multiqc_data
        │
        ├── trimmed/
        │   └── SRR30536566_full_DNA2_[R1,R2].trimmed.fastq.gz
        │
        ├── logs/
        │   └── cutadapt_SRR30536566_full_DNA2.log
        │   └── bwa_mem.log
        │   └── markduplicates.log
        │   └── SRR30536566_full_DNA2.flagstat.txt
        │   └── mutect2.stderr.log / mutect2.stdout.log
        │   └── learn_read_orientation_model.log
        │   └── get_pileup_summaries.log
        │   └── calculate_contamination.log
        │   └── filter_mutect_calls.log
        │   └── SRR30536566_full_DNA2.postfilter.log
        │
        ├── aligned/
        │   └── SRR30536566_full_DNA2.sorted.markdup.md.bam
        │   └── SRR30536566_full_DNA2.sorted.markdup.md.bam.bai
        │   └── SRR30536566_full_DNA2.markdup.metrics.txt
        │
        ├── variants/
        │   ├── *.contamination.table
        │   ├── *.f1r2.tar.gz
        │   ├── *.unfiltered.vcf.gz
        │   ├── *.unfiltered.vcf.gz.stats
        │   ├── *.unfiltered.vcf.gz.tbi
        │   ├── *.read-orientation-model.tar.gz
        │   ├── *.filtered.vcf.gz
        │   ├── *.filtered.vcf.gz.tbi
        │   ├── *.pileups.table
        │   ├── *.PASS.vcf.gz
        │   ├── *.PASS.vcf.gz.tbi
        │   └── *.postfiltered.vcf.gz
        │   └── *.postfiltered.vcf.gz.tbi
        │   └── *.postfilter_summary.txt
        │
        └── annotation/
```
---

## Results

### Comparison between outputs from `DNA` and `DNA2` environments using dataset `SRR30536566`

Pipeline outputs from the `DNA` and `DNA2` environments were compared using the log files generated during each pipeline step.

- Trimming and filtering: `cutadapt_SRR30536566_full_DNA2.log`

- alignment statistics: `bwa_mem.log`, `SRR30536566_full_DNA2.flagstat.txt`

- duplicate marking metrics `markduplicates.log`

- Mutect2 variant calling outputs: `mutect2.stderr.log`, `mutect2.stdout.log`

- LearnReadOrientationModel: `learn_read_orientation_model.log`

- GetPileupSummaries: `get_pileup_summaries.log`

- CalculateContamination: `calculate_contamination.log`

- FilterMutectCalls: `filter_mutect_calls.log`

- Variant postfilter: `SRR30536566_full.postfilter.log`

This comparison evaluates whether updating samtools in the `DNA2` environment has any measurable impact on alignment statistics, intermediate processing steps, and final somatic variant calls.


### Results - Table: Comparison of outputted metrics between files generated from environments `DNA` vs `DNA2`. Dataset: `SRR30536566`

| .log file                                    | Pipeline step            | Tool                      | Quantitative result  | DNA vs DNA2   |
| ---------------------------------------------| ------------------------ | ------------------------- | -------------------- | ------------- |
| `cutadapt_.log`                              | Trimming                 | Cutadapt                  | identical reads      | identical     |
| `bwa_mem.log`                                | Alignment                | BWA                       | 7,726,570 reads      | identical     |
| `bwa_mem.log`                                | Alignment runtime        | BWA                       | 881 s vs 711 s       | ~19% faster   |
| `markduplicates.log`                         | Duplicate marking        | Picard                    | 4,101,894 duplicates | identical duplicate detection |
| `SRR30536566_.flagstat.txt`                  | Mapping stats            | Samtools flagstat         | 99.32% mapped        | mapping identical;<br>minor differences due to<br>samtools reporting format     |
| `mutect2.stderr.log`<br>`mutect2.stdout.log` | Variant calling          | Mutect2                   | 948 candidate variants evaluated         | Identical variant-calling statistics     |
| `learn_read_orientation_model.log`           | Orientation bias model   | LearnReadOrientationModel | 32 sequence contexts modeled      | identical EM convergence and model     |
| `get_pileup_summaries.log`                   | Pileup summaries         | GetPileupSummaries        | 1,464,279 reads processed<br>901,738 filtered<br>88,955 loci analyzed          | identical     |
| `calculate_contamination.log`                | Contamination estimation | CalculateContamination    | 0 changepoints detected       | identical     |
| `filter_mutect_calls.log`                    | Variant filtering        | FilterMutectCalls         | 948 possible variant sites<br>237 variant candidates        | identical     |
| `SRR30536566_full.postfilter.log`            | Custom post-filter       | Panel thresholds          | 3 final variants     | identical     |

---

### Filtered VCF files comparison using **BCFtools**:

**I. Activate `DNA2` and go to `~/Genomics_cancer`**:

```bash
bcftools isec -p vcf_compare \
~/data/SRR30536566_full_DNA2/variants/SRR30536566_full_DNA2.filtered.vcf.gz \
~/data/SRR30536566_full/variants/SRR30536566_full.filtered.vcf.gz
```
> [!Note] Use absolute paths when using `bcftools isec -p vcf_compare`

Expected output: Generation of folder called `~/Genomics_cancer/vcf_compare`

```bash
-rw-r--r-- 1 Frano staff 181K Mar 10 08:34 0000.vcf
-rw-r--r-- 1 Frano staff 181K Mar 10 08:34 0001.vcf
-rw-r--r-- 1 Frano staff 261K Mar 10 08:34 0002.vcf
-rw-r--r-- 1 Frano staff 261K Mar 10 08:34 0003.vcf
-rw-r--r-- 1 Frano staff 5.7K Mar 10 08:34 sites.txt
-rw-r--r-- 1 Frano staff 1.5K Mar 10 08:36 README.txt
```

**Meaning**:

| File       | Meaning                                     |
| ---------- | ------------------------------------------- |
| `0000.vcf` | variants only in file 1 (**DNA2**)              |
| `0001.vcf` | variants only in file 2 (**DNA**)               |
| `0002.vcf` | variants from **DNA2** that are also present in **DNA**                     |
| `0003.vcf` | variants from **DNA** that are also present in **DNA2** |

**II. Count the amount of shared and intersected variants in `0000.vcf`, `0001.vcf`, `0002.vcf` and `0003.vcf`**

```bash
grep -vc "^#" 0000.vcf 
0
grep -vc "^#" 0001.vcf
0
grep -vc "^#" 0002.vcf
237
grep -vc "^#" 0003.vcf
237
```

**Meaning**:

| File       | Variants | Meaning                     |
| ---------- | -------- | --------------------------- |
| `0000.vcf` | **0**    | variants unique to **DNA2** |
| `0001.vcf` | **0**    | variants unique to **DNA**  |
| `0002.vcf` | **237**  | variants shared by both     |
| `0003.vcf` | **237**  | same shared variants        |

**Conclusion**: DNA variants = DNA2 variants. There are no variants unique to either pipeline, and that means `DNA` and `DNA2` environments produced **identical variant calls**.


**III. Check the file `sites.txt`**
`sites.txt` is a summary of variant sites produced by BCFtools.

| Column | Meaning                     |
| ------ | --------------------------- |
| 1      | Chromosome                  |
| 2      | Position                    |
| 3      | REF allele                  |
| 4      | ALT allele                  |
| 5      | The bitmask indicates in which<br>input file(s) the variant appears. |

Example:
```bash
chr1 114705278 A G 11
```
**Meaning**:

| Field     | Meaning                   |
| --------- | ------------------------- |
| chr1      | chromosome                |
| 114705278 | position                  |
| A         | reference                 |
| G         | alternative               |
| 11        | present in **both files** |

**Bitmask explanation**:

| Code | Meaning         |
| ---- | --------------- |
| `10` | only file 1 (**DNA2**)     |
| `01` | only file 2  (**DNA**)   |
| `11` | present in both |

**Conclusion**: The entire file shows: **11**, that means, for every variant → **all sites identical**.

### Visual summary

```bash
   DNA pipeline
        │
        ├── 237 variants
        │
        ▼
   DNA2 pipeline
        │
        ├── 237 variants
        │
        ▼
   Difference
        │
        └── 0 variants
```
**Alternative verification: In `~/Genomics_cancer`**: 

```bash
zgrep -v "^#" data/SRR30536566_full_DNA2/variants/SRR30536566_full_DNA2.filtered.vcf.gz | cut -f1,2 | sort | uniq | wc -l

zgrep -v "^#" data/SRR30536566_full/variants/SRR30536566_full.filtered.vcf.gz | cut -f1,2 | sort | uniq | wc -l
```
Ouput: 237 unique variants sites (in `DNA` and `DNA2`) → Difference = 0 (identical)

---

### Post-filtered VCF files comparison using **BCFtools**

Similarly, the comparison can be done with BCFtools from `*.postfiltered.vcf.gz`. Go to `~/Genomics_cancer`:

```bash
bcftools isec -p vcf_compare \
~/data/SRR30536566_full_DNA2/variants/SRR30536566_full_DNA2.postfiltered.vcf.gz \
~/data/SRR30536566_full/variants/SRR30536566_full.postfiltered.vcf.gz
```

Alternatively, go to `~/Genomics_cancer`

```bash
bcftools view -H data/SRR30536566_full/variants/SRR30536566_full.postfiltered.vcf.gz
bcftools view -H data/SRR30536566_full_DNA2/variants/SRR30536566_full_DNA2.postfiltered.vcf.gz
```

**Results - Table of comparison**: Both pipelines show these metrics
| Chr  | Pos       | REF | ALT | DP   | AD      | AF    | GT  | F1R2    | F2R1    | FAD     | SB              | POPAF | TLOD    |
| ---- | --------- | --- | --- | ---- | ------- | ----- | --- | ------- | ------- | ------- | --------------- | ----- | ------- |
| chr1 | 114713909 | G   | T   | 817  | 648,115 | 0.154 | 0/1 | 241,52  | 289,43  | 567,103 | 328,320,61,54   | 5.6   | 323.24  |
| chr3 | 179218294 | G   | A   | 1324 | 915,347 | 0.277 | 0/1 | 363,138 | 382,148 | 786,300 | 454,461,162,185 | 5.6   | 1026.18 |
| chr3 | 179226113 | C   | G   | 589  | 165,394 | 0.698 | 0/1 | 68,145  | 75,182  | 151,350 | 23,142,63,331   | 1.34  | 1387.91 |


> [!NOTE]
> - DP = total depth
> - AD = allelic depth (ref, alt)
> - AF = allelic fraction (ALT / DP)
> - GT = genotype (0/1 = heterozygous)
> - F1R2 / F2R1 – counts of alt alleles on forward/reverse strands (strand bias)
> - FAD – filtered allelic depth
> - SB – strand bias table (ref-forward, ref-reverse, alt-forward, alt-reverse)
> - POPAF – population allele frequency
> - TLOD – log-odds score for variant


**Comparison DNA vs DNA2**

Looking line by line:
- **chr1:114713909** G>T → identical

- **chr3:179218294** G>A → identical

- **chr3:179226113** C>G → identical

All fields, depth, AF, genotypes, filter tags **are exactly the same**.

> [!IMPORTANT]
> Mutect2 internally evaluates a large number of candidate variant events during its statistical modeling process.  
> In this dataset, **948 candidate variant records were evaluated**, but only **237 candidate variants were written to the VCF file** after applying the Mutect2 calling model.
>
> `FilterMutectCalls` then annotates these variants with filtering tags but does **not remove them from the VCF file**. This behavior is expected.
>
> When comparing the VCF files with **bcftools isec**, the comparison operates on explicit variant records defined by `CHROM + POS + REF + ALT`.
>
> The comparison confirms that **all 237 variants are identical between the `DNA` and `DNA2` environments**.
>
> Finally, a custom post-filtering step based on:
>
> - PASS status  
> - depth (DP)  
> - alternate allele count (AD)  
> - variant allele frequency (VAF)
>
> reduces the dataset to **3 high-confidence variants**.

```bash
 Mutect2
  │
  ├─ evaluates 948 possible sites
  │
  └─ writes 237 candidate variants → unfiltered.vcf
          │
          ▼
  FilterMutectCalls
          │
          └─ adds FILTER tags (still 237 variants)
                   │
                   ▼
      Custom pipeline thresholds
                   │
                   ├── Post-fitering
                   │  
                   ▼
      3 final variants in `DNA` and `DNA2`
```
---

### Quick way of checking VCF equality

**I. Go to `~/Genomics_cancer`**

```bash
zgrep -v "^#" data/SRR30536566_full/variants/SRR30536566_full.postfiltered.vcf.gz | md5sum
```
Expected output: `235266e81bb7ad44a73a1594cdd29291`

```bash
zgrep -v "^#" data/SRR30536566_full_DNA2/variants/SRR30536566_full_DNA2.postfiltered.vcf.gz | md5sum
```
Expected output: `235266e81bb7ad44a73a1594cdd29291`

**Meaning**: If the hashes match, the variant records are identical. This is **much faster** than `bcftools isec`, especially for large VCFs. Header lines are excluded using `zgrep -v "^#"`.

**Alternative**

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' data/SRR30536566_full_DNA2/variants/SRR30536566_full_DNA2.postfiltered.vcf.gz | md5sum
```
Expected output: `f085e93f16b007dd2f3b0a13875b405b`

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' data/SRR30536566_full/variants/SRR30536566_full.postfiltered.vcf.gz | md5sum
```
Expected output: `f085e93f16b007dd2f3b0a13875b405b`

---

## Conclusion

The comparison between the `DNA` and `DNA2` environments shows that the final post-filtered VCF files are **identical**.

No differences were observed in:

- variant positions
- reference and alternate alleles
- genotypes
- depth and allele frequency metrics
- filtering annotations

This confirms that upgrading **samtools** from version **0.1.19** to **1.22.1** does **not affect the somatic variant calls produced by the pipeline**.

Across all intermediate processing steps and final outputs, the results generated by `DNA` and `DNA2` are fully consistent.

Therefore, the updated `DNA2` environment can safely replace the legacy `DNA` environment while maintaining full reproducibility of the analysis.

---

Top of page 👉 [Part V: Pipeline maintenance and Environment Validation](README_Part5_DNA2_pipeline_update.md##part-v--pipeline-maintenance-and-environment-validation)

Previous analysis in `DNA` environment (with legacy samtools) 👉 [Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline](README_Part4_fullbash.md)

Jump to the first part of this tutorial 👉 [Part I – Preparation & setup](README_setup_Part1-3.md)

Main page 👉 [Bash_pipeline_NGS](README.md)

