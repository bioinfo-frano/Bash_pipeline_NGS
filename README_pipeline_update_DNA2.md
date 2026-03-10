# Part V – Pipeline maintenance

## Table of Contents

- [Introduction](#introduction)
- [Methodology](#methodology)
    - [Environment reproducibility (`DNA2`)](#3-environment-reproducibility-dna2)
    - [Preparation of folder structure and update of Bash pipeline](#5-preparation-of-folder-structure-and-update-of-bash-pipeline)
- [Results](#results)
    - [Verify the differences between the filtered VCF files using **BCFtools**](#verify-the-differences-between-the-filtered-vcf-files-using-bcftools)
    - [Verify the differences between the post-filtered VCF files](#verify-the-differences-between-the-post-filtered-vcf-files)
- [Conclusion](#conclusion)

# Introduction

Maintaining or updating a bioinformatics pipeline is a critical responsibility for any bioinformatician, especially when working with research or clinical genomic data.

Pipeline maintenance ensures that the code written months or years ago still produces **reliable, reproducible, and scientifically valid results today**.

Because bioinformatics software evolves constantly, pipelines can gradually become outdated. New tool versions, new reference datasets, and improvements in algorithms may affect analysis results. Without proper maintenance, a pipeline may eventually stop working or produce inconsistent results.

Therefore, maintaining a pipeline throughout its lifecycle is an essential part of bioinformatics practice.

In general, pipeline maintenance involves several types of tasks.

## 1. Tool and Dependency updates

Bioinformatics software evolves rapidly. Tools such as GATK, BWA and samtools are regularly updated, along with their dependencies (e.g. Python, R, or Java). For example, if a pipeline originally used **GATK 4.2** and is later updated to **GATK 4.5**, it is important to verify that the updated pipeline still produces consistent results.

Another important aspect of maintenance is ensuring **reproducibility of the computational environment**. Package managers such as **Conda** allow the exact software environment used in the past to be recreated later. In more advanced pipelines, container technologies such as Docker or Singularity may also be used to fully reproduce the runtime environment.

## 2. Exceptional cases

Occasionally, a specific sequencing sample may produce an unexpected error during pipeline execution.

A robust pipeline should ideally detect these failures and allow the rest of the samples to continue processing. In practice, this may require adjusting the script to handle exceptional cases without interrupting the entire analysis run.

## 3. Reference Data updates

Many genomic reference resources are continuously updated, including:

- dbSNP

- ClinVar

- gnomAD

- Human reference genome GRCh38.

Maintaining a pipeline may therefore require downloading updated reference datasets, indexing them, and validating that they are compatible with the pipeline before sharing them with collaborators.

## 4. Bug fixes

During routine pipeline usage, unexpected issues may arise such as:

- incorrect file paths

- software incompatibilities

- syntax changes in tool commands

- incorrect parameter usage

Fixing these issues requires updating the pipeline code and validating that the modifications do not alter the expected biological results.

## 5. Pipeline Versioning and Provenance Tracking

Another essential aspect of pipeline maintenance is **tracking which version of the pipeline produced each result**. This concept is known as **pipeline versioning** and **data provenance**.

Even small changes in a pipeline may affect results. Examples include:

- updating samtools

- updating the variant caller Mutect2 from the GATK toolkit

- modifying filtering thresholds

- updating annotation databases such as gnomAD or ClinVar

To track these changes, bioinformatics pipelines typically use **version control systems** such as:

- Git

Code repositories such as GitHub allow developers to track every modification made to the pipeline.

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

# Practical example of pipeline maintenence

In this tutorial, the somatic variant analysis from a targeted NGS gene panel is performed using the Bash pipeline:

[`09_full_somatic_NGS_bash_script.sh`](bash_scripts/09_full_somatic_NGS_bash_script.sh). 

The script sequentially runs several bioinformatics tools, including the variant caller **Mutect2** from the GATK toolkit. All required software is installed inside a Conda environment called `DNA`. One of the tools used in this environment is **samtools**, which is widely used for manipulating sequencing alignment files (BAM/CRAM).

The original `DNA` environment contains a very old version of samtools (**0.1.19**, a "legacy" version). Modern versions of samtools (>1.10) rely on the library HTSlib and include numerous improvements and bug fixes.

An interesting question is therefore:

"***Will the number or identity of detected variants change if the pipeline is executed using a modern version of samtools instead of the legacy version?***"

---

## Methodology

To evaluate the potential impact of this update, the following approach will be used.

### 1. Compare the environments `DNA` and `DNA2`

A new environment called `DNA2` will be created with a modern version of samtools.

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

In this experiment, the primary tool intentionally modified is **samtools**, which is updated from version **0.1.19** to **1.22.1**. Because modern samtools releases are built on the **HTSlib** library, updating samtools also introduces compatible versions of HTSlib and related tools such as bcftools. 

All other core tools remain unchanged to isolate the effect of this update.

---

### 2. Create the new conda environment `DNA2`. 

The installation of `DNA2` should be performed from the `base` environment. All Conda environments are created from the `base` installation but remain isolated from each other.

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

This configuration ensures a **strict channel priority**, which helps prevent dependency conflicts when installing bioinformatics software.  

> [!CAUTION]
> During environment creation, dependency conflicts may occur when incompatible versions of tools are requested simultaneously. In such cases, removing strict version constraints often allows Conda to automatically resolve compatible versions. In genomics pipelines this frequently occurs with legacy Perl-based tools such as **Ensembl-VEP**. For this reason, complex pipelines often split tools into multiple environments to avoid dependency conflicts.

**III. Create the environment**:

```bash
conda create -n DNA2 \c 
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

This environment will likely produce dependency conflicts during installation. The conflict occurs because **Ensembl-VEP depends on the legacy samtools API (<0.2)** through `perl-bio-samtools`, preventing Conda from installing modern samtools versions (>=1.x) within the same environment.

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

**Alternative**: Instead of Conda, the environment can be created using mamba, which resolves dependencies much faster.

Check if mamba is installed:

```bash
which mamba
mamba --version
```

Then create the environment:

```bash
`mamba create -n DNA2 ...
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

**I. Export the environment**

After creating `DNA2`

```bash
conda activate DNA2
conda env export --no-builds > DNA2_conda_environment.yml
```

The file `DNA2_conda_environment.yml` will store all packages and versions of `DNA2`.

To get a higher level of reproducibility, remove `--no-builds` from the code.


**II. Create a fully reproducible lock file**

If you want a perfect reproducibility, then the strongest methods are:

```bash
conda list --explicit > DNA2_conda_environment.lock
```

or 

```bash
# First export DNA2 env without builds
conda env export --no-builds > DNA2_conda_environment.yml
# Second, create the lock DNA2 env specifically for macOS(Intel)
conda-lock lock -f DNA2_conda_environment.yml -p osx-64
# Alternatively, for macOS(Apple Silicon)
conda-lock lock -f DNA2_conda_environment.yml -p osx-arm64
# Or if you also want to increase portability to other platforms (linux-64, win-64)
conda-lock -f DNA2_conda_environment.yml
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


**III. Recreate the environment anywhere from the lock file**

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

### 4. Verifying how "clean" is `DNA2` environment

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
**Meaning** of `All requested packages already installed`: It means that samtools, htslib, bcftools were already installed in `DNA2` environment, and there is **no dependency resolution needed**.


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

### 6. Run the pipeline in `DNA2` environment

**I. Activate the environment**

```bash
conda activate DNA2
```

**II. Move to ~/scripts and run the `09_full_somatic_DNA2_updated.sh`**

```bash
cd ~/Genomics_cancer/scripts
./09_full_somatic_DNA2_updated.sh
```
The entire pipeline will the same sequencing dataset `SRR30536566`.

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

The resulting quality metrics will be compared using the reports recorded in `.logs`, paying particular attention to:

- Trimming and filtering: `cutadapt_SRR30536566_full_DNA2.log`

- alignment statistics: `bwa_mem.log`, `SRR30536566_full_DNA2.flagstat.txt`

- duplicate marking metrics `markduplicates.log`

- Mutect2 variant calling outputs: `mutect2.stderr.log`, `mutect2.stdout.log`

- LearnReadOrientationModel: `learn_read_orientation_model.log`

- GetPileupSummaries: `get_pileup_summaries.log`

- CalculateContamination: `calculate_contamination.log`

- FilterMutectCalls: `filter_mutect_calls.log`

- Variant postfilter: `SRR30536566_full.postfilter.log`

This comparison allows us to evaluate whether updating samtools has any measurable impact on the final somatic variant calls.


### Results - Table: Comparison of outputted metrics between files generated from environments `DNA` vs `DNA2`. Dataset: `SRR30536566`

| .log file                                    | Pipeline step            | Tool                      | Quantitative result  | DNA vs DNA2   |
| ---------------------------------------------| ------------------------ | ------------------------- | -------------------- | ------------- |
| `cutadapt_.log`                              | Trimming                 | Cutadapt                  | identical reads      | identical     |
| `bwa_mem.log`                                | Alignment                | BWA                       | 7,726,570 reads      | identical     |
| `bwa_mem.log`                                | Alignment runtime        | BWA                       | 881 s vs 711 s       | ~19% faster   |
| `markduplicates.log`                         | Duplicate marking        | Picard                    | 4,101,894 duplicates | identical duplicate detection |
| `SRR30536566_.flagstat.txt`                  | Mapping stats            | Samtools flagstat         | 99.32% mapped        | mapping identical;<br>minor differences due to<br>samtools reporting format     |
| `mutect2.stderr.log`<br>`mutect2.stdout.log` | Variant calling          | Mutect2                   | 948 variants         | Identical variant-calling statistics     |
| `learn_read_orientation_model.log`           | Orientation bias model   | LearnReadOrientationModel | 32 sequence contexts modeled      | identical EM convergence and model     |
| `get_pileup_summaries.log`                   | Pileup summaries         | GetPileupSummaries        | 1,464,279 reads processed<br>901,738 filtered<br>88,955 loci analyzed          | identical     |
| `calculate_contamination.log`                | Contamination estimation | CalculateContamination    | 0 changepoints detected       | identical     |
| `filter_mutect_calls.log`                    | Variant filtering        | FilterMutectCalls         | 948 variants retained         | identical     |
| `SRR30536566_full.postfilter.log`            | Custom post-filter       | panel thresholds          | 3 final variants     | identical     |

---

### Verify the differences between the filtered VCF files using **BCFtools**:

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
| `0000.vcf` | variants only in file 1 (DNA2)              |
| `0001.vcf` | variants only in file 2 (DNA)               |
| `0002.vcf` | variants shared by both                     |
| `0003.vcf` | variants present in both (intersection set) |

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
| 5      | bitmask indicating presence |

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
| `10` | only file 1     |
| `01` | only file 2     |
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

---

### Verify the differences between the post-filtered VCF files

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

**Meaning**: If the hashes match, the files are similar. This is **much faster** than `bcftools isec`, especially for large VCFs.

---

## Conclusion:

After post-filtering, the final VCFs from `DNA` and `DNA2` are identical.

- There is no difference in variants, genotypes, or quality metrics.

- This fully confirms that switching to `DNA2` is reproducible at the final variant call level.

- By comparing the generated metrics from the somatic analysis done in both Conda environments, the new environment **`DNA2` is a fully reproducible replacement of `DNA`**.

---

Go back to the top of 👉 [Part V: Pipeline maintenance](README_pipeline_update_DNA2.md#part-v--pipeline-maintenance)

Jump to the first part of this tutorial 👉 [Part I – Preparation & setup](README_setup_Part1-3.md)

Go to the main page 👉 [Bash_pipeline_NGS](README.md)

