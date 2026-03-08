# Pipeline maintenance

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

# Practical example

In this tutorial, the somatic variant analysis from a targeted NGS gene panel is performed using the Bash pipeline:

[`09_full_somatic_NGS_bash_script.sh`](bash_scripts/09_full_somatic_NGS_bash_script.sh). 

The script sequentially runs several bioinformatics tools, including the variant caller **Mutect2** from the GATK toolkit. All required software is installed inside a Conda environment called `DNA`. One of the tools used in this environment is **samtools**, which is widely used for manipulating sequencing alignment files (BAM/CRAM).

The original `DNA` environment contains a very old version of samtools (**0.1.19**, a "legacy" version). Modern versions of samtools (>1.10) rely on the library HTSlib and include numerous improvements and bug fixes.

An interesting question is therefore:

"***Will the number or identity of detected variants change if the pipeline is executed using a modern version of samtools instead of the legacy version?***"

---

## Methodology: 

To evaluate the potential impact of this update, the following approach will be used.

### 1. Compare the environments `DNA` and `DNA2`

A new environment called `DNA2` will be created with a modern version of samtools and updated dependencies.

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

In this experiment, the primary tool intentionally modified is samtools, which is updated from version 0.1.19 to 1.22.1. Because modern samtools releases are built on the HTSlib library, updating samtools also introduces compatible versions of HTSlib and related tools such as bcftools. All other core tools remain unchanged to isolate the effect of this update.

---

### 2. Create the new conda environment `DNA2`. 

The installation of `DNA2` should be in `base` environment. All Conda environments are created from the `base` installation but remain isolated from each other.

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

You should see something like:

```bash
channels:
  - conda-forge
  - bioconda
  - defaults
```

This provides a strict channel priority that helps prevent dependency conflicts when installing bioinformatics software.  

> WARNING:
> During environment creation, dependency conflicts may occur when incompatible versions of tools are requested simultaneously. In such cases, removing strict version constraints often allows Conda to automatically resolve compatible versions. In genomics pipelines this frequently occurs with legacy Perl-based tools such as **Ensembl-VEP**. For this reason, complex pipelines often split tools into multiple environments to avoid dependency conflicts.

**III. Then run**:

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

This environment will likely produce dependency conflicts during installation. The conflict is basically, as it was mentioned, due to the dependency `ensembl-vep`.
`Ensembl-VEP` depends on the Perl module `perl-bio-samtools`, which relies on the legacy samtools API (<0.2), preventing Conda from installing modern samtools versions (>=1.x)
within the same environment.

```bash
ensembl-vep
  └─ perl-bioperl
       └─ perl-bio-samtools
            └─ samtools (<0.2 legacy API)
```

However, we need `samtools 1.22.1`. To resolve this conflict, one option is to remove `ensembl-vep` and let Conda decide, which `Perl` version should be installed automatically.

The installation of `DNA2` can therefore be repeated with the following modifications.

IV. Recreate the environment without Ensembl-VEP:

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

It will be shown in terminal:

```bash
Solving environment: done
...
Executing transaction: done
```

>**Key-Note**: Alternatively, instead of creating `DNA2` with Conda, it can be use a better (faster) alternative with **mamba**.
In base environment, check whether **mamba** is installed with:
>
> `which mamba`
>
> `mamba --version`
>
> In order to create `DNA2` use the same Conda command replacing `conda` for `mamba` like this:
>
> `mamba create -n DNA2 ...`
>
> Then activate:
>
> `conda activate DNA2`
>
> `conda list`  

So now, the new `DNA2` environment will have a newer version of samtools.

| Dependency | Version (`DNA`) | Dependency   | Version (`DNA2`) |
| ---------- | --------------- | ------------ | ---------------- |
| samtools   | 0.1.19          | **samtools** | **1.22.1**       |

---

### 3. Environment reproducibility (`DNA2`)

Instead of installing environments manually, you can export a full environment using a `.yml` file.

**I. Export the environment**

After creating `DNA2`

```bash
conda activate DNA2
```

Run:

```bash
conda env export --no-builds > DNA2_conda_environment.yml
```

This will create `DNA2_conda_environment.yml`


**II. Recreate the environment anywhere**

Anyone can recreate the **exact same pipeline environment**:

```bash
conda env create -f DNA2_conda_environment.yml
```

If the `.yml` should go directly to a folder `/envs` then:

```bash
conda env create -f envs/DNA2_conda_environment.yml
```

**III. Alternative method**: Lock your environment.

```bash
conda list --explicit > DNA2_lock.txt
```

Then recreate byte-identical environments:

```bash
conda create --name DNA2 --file DNA2_lock.txt
```

This is **maximum reproducibility**.


These are **standard methods used in bioinformatics to transfer exact conda environments**. With the **YAML** file, the environment won't be affected by changes in versions on each dependency, making the environment:

- reproducible
- shareable
- version-controlled

>IMPORTANT:
>
> By creating the environment from a **YAML** file, i.e. `conda env create -f envs/DNA2_conda_environment.yml`, Conda will **resolve dependencies again** and install compatible packages that satisfy the constraints. This means that when installing the environment through `.yml` in another system, Conda might install different builds. For example, `zlib 1.3.1` vs `zlib 1.3.1 build_1`. So the environment is reproducible **conceptually**, but not identical.
> On the other hand, by creating the environment from a **lock** file, i.e. `conda create --name DNA2 --file DNA2_lock.txt`, Conda uses an **explicit package list**. In this case, Conda **does not solve dependencies**, and it installs the **exact package builds**. So you get **bit-identical environments**. For example, `samtools-1.22.1-h96c455f_0` will be the same version, same build, same hash. The problem with this way of transferring environment is that they are **platform specific**, for example **Platform: osx-64**. This means that the installation of the environment **would fail** on:
> - Linux cluster
> - Apple Silicon
> - HPC systems
>
> YAML files **are portable**, **reproducible**, and **standard in bioinformatics**.

---

### 4. Tools verification 

After activating the environment, it is recommended to verify that the main tools are correctly installed:

```bash
samtools --version
bcftools --version
gatk --version
```

---

### 5. Preparation of folder structure and Update of Bash pipeline

a) Copy the bash script `09_full_somatic_NGS_bash_script.sh` and create a `09_full_somatic_DNA2_updated.sh`

Go to ~/Genomics_cancer/scripts

```bash
cp 09_full_somatic_NGS_bash_script.sh 09_full_somatic_DNA2_updated.sh
```

b) Create a new folder called "SRR30536566_full_DNA2" with its subfolders:

Go to ~/Genomics_cancer/data

```bash
mkdir -p SRR30536566_full_DNA2/{aligned,logs,qc,trimmed,variants,annotation}
```

Folder structure:

```code
Genomics_cancer/
├── data/
│   ├── SRR30536566_full/   # Used by the unified script. DNA environment.
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── logs/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   └── SRR30536566/        # Used for step-by-step bash scritp workflow (Part I & II)
│   │   ├── raw_fastq/
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   └── SRR30536566_full_DNA2/        # Used by the unified script. DNA2 environment.
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

c) Update the bash script `09_full_somatic_DNA2_updated.sh`

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

### 6. Run the pipeline in the DNA2 environment

Activate the environment:

```bash
conda activate DNA2
```

Move to scripts (where `09_full_somatic_DNA2_updated.sh` is)

```bash
cd ~/Genomics_cancer/scripts
```

Run:

```bash
bash 09_full_somatic_DNA2_updated.sh
```

or 

```bash
./ 09_full_somatic_DNA2_updated.sh
```

### Expected output structure after the run

```code
Genomics_cancer/
└── data/
    └── SRR30536566_full_DNA2/
        ├── qc/
        │   ├── raw/
        │   ├── trimmed/
        │   └── md_flagstat/
        │
        ├── trimmed/
        │
        ├── logs/
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

The entire pipeline will then be executed again on the same sequencing dataset (`SRR30536566`)

The resulting quality metrics will be compared using the report generated by MultiQC, paying particular attention to:

- alignment statistics

- duplicate marking metrics

- the number of variants detected

- the number of variants remaining after filtering

This comparison allows us to evaluate whether updating samtools has any measurable impact on the final somatic variant calls.

MAKE A TABLE COMPARING OUTPUTS from OLD `DNA` vs NEW `DNA2`


