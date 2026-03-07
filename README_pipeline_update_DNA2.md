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

### 3. Reference Data updates

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


# Practical example

In this tutorial, the somatic variant analysis from a targeted NGS gene panel is performed using the Bash pipeline (`09_full_somatic_NGS_bash_script.sh` (make a link)). The script sequentially runs several bioinformatics tools, including the variant caller Mutect2 from the GATK toolkit.

All required software is installed inside a Conda environment called `DNA`. One of the tools used in this environment is samtools, which is widely used for manipulating sequencing alignment files (BAM/CRAM).

The original `DNA` environment contains a very old version of samtools (0.1.19, "legacy"). Modern versions of samtools (>1.10) rely on the library HTSlib and include numerous improvements and bug fixes.

An interesting question is therefore:

"***Will the number or identity of detected variants change if the pipeline is executed using a modern version of samtools instead of the legacy version?***"


## Methodology: 

To evaluate the potential impact of this update, the following approach will be used.

### 1. Compare "DNA" with a new Conda environment.

A second environment called `DNA2` will be created. This environment will contain a more recent version of samtools and updated dependencies (for example **Perl 5.32.1**).


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

MAKE A TABLE CALLED: "CONDA ENV 'OLD' (`DNA`) and NEW (`DNA2`) COMPARISON". IN THIS TABLE, THE FIRST TWO COLUMNS SHOWS THE DEPENDENCIES FROM OLD AND VERSION, AND THE LAST TWO THE DEPENDENCIES OF THE NEW AND VERSION

### 2. Create the new conda environment `DNA2`. 


```bash
conda create -n DNA2 \
  -c conda-forge -c bioconda -c defaults \
  python=3.9 \
  openjdk=17 \
  perl=5.32 \
  fastqc \
  multiqc \
  cutadapt \
  bwa \
  samtools***VERSION \
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

CORRECT THE CODE



### 3. Update the Bash pipeline

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

c) Update the bash script `09_full_somatic_DNA2_updated.sh`

- Change variables in "Configuration"

SAMPLE="SRR30536566"

Old version
`SAMPLE_FULL="${SAMPLE}_full"`

New version
`SAMPLE_FULL="${SAMPLE}_full_DNA2"`

- Change samtools syntax in "Step3: Sort BAM (coordinate sort)"

Old alignment block
```bash
samtools sort -@ "$THREADS" \
  "$ALIGN_DIR/${SAMPLE_FULL}.bam" \
  "$ALIGN_DIR/${SAMPLE_FULL}.sorted"
```

New version
```bash
samtools sort -@ "$THREADS" \
  "$ALIGN_DIR/${SAMPLE_FULL}.sorted" \
  "$ALIGN_DIR/${SAMPLE_FULL}.bam"
```


### 3. Re-run the pipeline

The entire pipeline will then be executed again on the same sequencing dataset (`SRR30536566`)

The resulting quality metrics will be compared using the report generated by MultiQC, paying particular attention to:

- alignment statistics

- duplicate marking metrics

- the number of variants detected

- the number of variants remaining after filtering

This comparison allows us to evaluate whether updating samtools has any measurable impact on the final somatic variant calls.

MAKE A TABLE COMPARING OUTPUTS from OLD `DNA` vs NEW `DNA2`


