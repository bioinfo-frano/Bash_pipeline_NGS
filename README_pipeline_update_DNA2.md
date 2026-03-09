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

The installation of `DNA2` should be performed from the in `base` environment. All Conda environments are created from the `base` installation but remain isolated from each other.

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

> [!WARNING]
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
Expected output: `DNA2_conda_environment.yml` in workig directory

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
conda create --name DNA2 --file env_lock.txt
```
or
```bash
conda create --n DNA2 --file env_lock.loc
```

These are **standard methods used in bioinformatics to transfer exact Conda environments**.

> [!IMPORTANT]
> By creating the environment from a **YAML** file, i.e. `conda env create -f envs/DNA2_conda_environment.yml`, Conda will **resolve dependencies again** and install compatible packages that satisfy the constraints. This means that when installing the environment through `.yml` in another system, Conda might install different builds. For example, `zlib 1.3.1` vs `zlib 1.3.1 build_1`. So the environment is reproducible **conceptually**, but not identical.
> On the other hand, by creating the environment from a **lock** file, i.e. `conda create --name DNA2 --file DNA2_lock.txt`, Conda uses an **explicit package list**. In this case, Conda **does not solve dependencies**, and it installs the **exact package builds**. So you get **bit-identical environments**. For example, `samtools-1.22.1-h96c455f_0` will be the same version, same build, same hash. The problem with this way of transferring environment is that they are **platform specific**, for example **Platform: osx-64**. This means that the installation of the environment **would fail** on:
> - Linux cluster
> - Apple Silicon
> - HPC systems
>
> YAML files **are portable**, **reproducible**, **shareable**,and **standard in bioinformatics**.


### SUMMARY

`conda list --explicit > DNA2.lock`  ➡️ **Creates a lock file from an existing environment**. 

What it does:

- Exports exact package URLs

- Includes build numbers

- Includes channels

- Includes platform

Example inside the file:

```bash
@EXPLICIT
https://conda.anaconda.org/bioconda/osx-64/samtools-1.22.1-h96c455f_0.conda
https://conda.anaconda.org/bioconda/osx-64/bcftools-1.22-hb1a7a94_0.conda
```
Install these **exact binaries**, no solving. **This is the most reproducible export**.


`conda create -n DNA2 --file DNA2.lock`  ➡️ **Recreates the exact environment and save it to a file**. 

Important:

- No dependency solving

- Installs exact builds

- Same packages

- Same versions

- Same builds

This guarantees **bit-identical environments**. Used for HPC pipelies, published analyses and reproducibility


`conda env create -f DNA2_conda_environment.yml` or `conda create --name DNA2_clone --file env.txt`: both do not have option `--explicit`

Then Conda must:

- solve dependencies again

- choose builds again

- maybe upgrade things

So the environment may be **different**.

**IMPORTANT**: There's no technical difference in filenames:

```bash
DNA2.lock
DNA2_lock.txt
env.txt
```
What matters is the **content**, especially the line:
```bash
@EXPLICIT
```

`conda-lock -f DNA2_conda_environment.yml` ➡️ **These contain fully solved environments**.

---

### 4. Tools verification 

After activating the environment, it is recommended to verify that the main tools are correctly installed:

```bash
samtools --version
bcftools --version
gatk --version
```

---

### 5. Verifying how "clean" is `DNA2` environment

A “clean” environment means that all packages are resolved correctly with no hidden conflicts.

Before using the environment `DNA2`:

**1. Activate the new environment**:

```bash
conda activate DNA2
```

**2. Check for broke dependencies**:

```bash
conda list --explicit
```
output:

```bash
# This file may be used to create an environment using:
# $ conda create --name <env> --file <this file>
# platform: osx-64
@EXPLICIT
https://conda.anaconda.org/conda-forge/osx-64/coreutils-9.5-h10d778d_0.conda
https://conda.anaconda.org/conda-forge/noarch/_r-mutex-1.0.1-anacondar_1.tar.bz2
...
https://conda.anaconda.org/conda-forge/noarch/seaborn-base-0.13.2-pyhd8ed1ab_3.conda
https://conda.anaconda.org/bioconda/noarch/picard-3.4.0-hdfd78af_0.tar.bz2
https://conda.anaconda.org/conda-forge/noarch/seaborn-0.13.2-hd8ed1ab_3.conda
```
**Meaning**: This codes is useful when creating a **lock file** listing the exact packages. This means that it's possible to recreate the **identical environment** later using:

```bash
conda create --name DNA2_clone --file env.txt
```
where `env.txt` is the saved file. This guarantees **bit-identical environments**.

Alternatively, **freeze the environment** like this:

```bash
conda list --explicit > DNA2.lock
```
Then the environment can be recreated exactly on another machine.

Example:

```bash
conda create -n DNA2 --file DNA2.lock
```
This is important when installing an environment in an HPC cluster, reproducing papers and pipeline sharing.


**3. Check the installation history**:

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



**4. Simulate installing critical packages to check for conflicts**:

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


>**Key note**: Using `--dry-run` ensures that you detect dependency conflicts **without modifying the environment**, preventing surprises later in your analysis.


---

### 6. Preparation of folder structure and Update of Bash pipeline

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


