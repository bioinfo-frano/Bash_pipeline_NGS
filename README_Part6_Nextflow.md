# Part IV – Nextflow Pipeline: Fully Automated Somatic DNA-NGS Pipeline (Single Nextflow Script)

In [Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline](README_Part4_fullbash.md) it was introduced a single bash script with the full NGS somatic analysis pipeline of the 7-gene amplicon dataset `SRR30536566`. This time is the turn of showing a Nextflow pipeline script of the same dataset, which can also do the complete analysis of the same dataset.

## What is Nextflow?

**Nextflow** is an open-source, data-driven workflow management system designed for creating scalable, portable, and reproducible computational pipelines, primarily in bioinformatics. The scripting language is not really trivial, and it takes a while to get know with its uncommon concepts, e.g. nextflow is written in groovy. 
Nextflow can deploy workflows on a variety of execution platforms, including your local machine, HPC schedulers, and cloud.

## What are the advantages and disadvantages of Nextflow and Bash pipelines?

| **Aspect**                | **Bash Script**                                                                                                                                                                                                 | **Nextflow Pipeline**                                                                                                                                                                                                                     |
|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Scalability**           | ❌ **Limited to one sample** – designed for a single sample; running multiple samples requires manual parallelisation or wrapper scripts.                                                                         | ✅ **Built‑in multi‑sample support** – channel logic easily handles any number of samples; processes run in parallel automatically.                                                                                                       |
| **Reproducibility**       | ⚠️ **Moderate** – relies on system‑wide tools and manual version management; environment may drift over time.                                                                                                    | ✅ **High** – uses Conda environment (`DNA2`) with exact tool versions; environment is self‑contained and reproducible across machines.                                                                                                    |
| **Portability**           | ❌ **Low** – paths are hard‑coded; moving to another system requires editing many lines.                                                                                                                         | ✅ **High** – all paths are defined as parameters at the top; changing them once adapts the whole pipeline.                                                                                                                               |
| **Error handling & resumption** | ❌ **Manual** – if a step fails, the whole script stops; you must fix and restart from the beginning (unless you manually add checkpoints).                                                              | ✅ **Automatic** – each step is isolated; after a fix, `-resume` continues from the last successful process, saving time and compute.                                                                                                     |
| **Logging & output organization** | ⚠️ **Ad‑hoc** – log files are written to specific directories, but organisation depends on manual `mkdir` and naming conventions.                                                                             | ✅ **Structured** – `publishDir` directives place outputs in clear, consistent folders (`qc/`, `logs/`, `aligned/`, `variants/`). Logs are automatically collected.                                                                       |
| **Resource management**   | ❌ **None** – all steps use the same thread/memory settings; no way to request different resources per task.                                                                                                    | ✅ **Per‑process control** – each process can specify `cpus`, `memory`, and `time`; resources are allocated appropriately (e.g., BWA‑MEM gets fewer threads to reduce memory).                                                            |
| **Parallel execution**    | ❌ **Sequential** – steps run one after another, even if they could run in parallel (e.g., FastQC on raw and trimmed reads).                                                                                    | ✅ **Parallel** – independent processes run concurrently, reducing total runtime (e.g., trimming and alignment run in parallel with QC steps).                                                                                            |
| **Code modularity & reuse** | ❌ **Monolithic** – a single long script; modifying a step (e.g., changing alignment parameters) risks breaking other parts.                                                                                  | ✅ **Modular** – each step is a separate process; you can update, replace, or add processes without affecting others.                                                                                                                     |
| **Dependency management** | ❌ **Manual** – you must ensure all tools (`bwa`, `samtools`, `gatk`, etc.) are installed and in `PATH`; conflicts can arise.                                                                                  | ✅ **Automatic** – the Conda environment (`DNA2`) bundles all tools with correct versions; no external dependencies beyond Conda.                                                                                                         |
| **Learning curve**        | ✅ **Low** – simple, linear script; easy for anyone with basic Bash knowledge to understand and modify.                                                                                                          | ⚠️ **Steeper** – requires understanding of Nextflow concepts (channels, processes, DSL2, etc.). However, once learned, it is much more powerful.                                                                                        |
| **Debugging**             | ⚠️ **Straightforward but manual** – you can run individual commands from the script, but tracking down errors in a long script can be tedious.                                                                   | ✅ **Isolated** – each process runs in its own work directory with `.command.sh`, `.command.log`, etc.; you can inspect exactly what happened for a single step.                                                                          |
| **Workflow visibility**   | ❌ **Low** – the overall flow is implicit; you have to read the whole script to understand the pipeline.                                                                                                        | ✅ **High** – the `workflow` block clearly shows the data flow between processes, making the pipeline self‑documenting.                                                                                                                  |
| **Job scheduling / HPC readiness** | ❌ **Not applicable** – runs only on a local machine; no integration with cluster schedulers (SLURM, PBS, etc.).                                                                                          | ✅ **Flexible** – Nextflow supports local execution, but can also submit jobs to HPC clusters or cloud with minimal configuration changes (`-with-sge`, `-with-slurm`, etc.).                                                            |
| **Development time**      | ✅ **Quick to write** – a linear script can be written and tested rapidly for a single sample.                                                                                                                   | ⚠️ **Longer initial development** – requires designing channels and processes, but pays off with reusability and robustness.                                                                                                             |
| **Final output consistency** | ✅ **Matches expectations** – produces exactly the files and logs you need, as we verified.                                                                                                                      | ✅ **Matches exactly** – the Nextflow pipeline was designed to mirror the Bash script, so outputs are identical (and even better organised).                                                                                              |


### References: 

- **Nextflow documentation**: <https://www.nextflow.io/docs/latest/index.html>

- **Maxwell Cluster**: <https://wiki.desy.de/maxwell/documentation/workflows/nextflow/>

- **nf-core**: "*A global community collaborating to build open-source Nextflow components and pipelines*" <https://nf-co.re/>


## Nextflow analysis of dataset `SRR30536566`

## 1. Install Nextflow to conda environment 'DNA'

**In Terminal**

1. Go to DNA environment

```bash
conda activate DNA
```

⚠️ IMPORTANT: Check whether DNA environment has already Java version >= 17 installed. 

In Terminal:

```bash
conda list
```

Check the version of Java manually by scrolling and finding the dependency "openjdk"

Alternative:

```bash
java -version
```
Output:

```bash
java -version
openjdk version "17.0.17" 2025-10-21 LTS
```

2. Nextflow installation (official way)

```bash
curl -s https://get.nextflow.io | bash
```

Output:

```bash
N E X T F L O W
      version 25.10.4 build 11173
      created 10-02-2026 15:17 UTC (16:17 CEST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io

Nextflow installation completed. Please note:
- the executable file `nextflow` has been created in the folder: /Users/Frano/Desktop/Bioinfo_2026/Genomics_cancer/scripts
- you may complete the installation by moving it to a directory in your $PATH
```

***This command is not installing Nextflow via conda, just downloading an executable script called "nextflow" in your current directory. That means that by typing "conda list", you won't find "nextflow"***
 

3. Now move it into your conda environment’s bin:

```bash
mv nextflow $CONDA_PREFIX/bin/
chmod +x $CONDA_PREFIX/bin/nextflow
```

4. Verify installation:

```bash
nextflow -version
```
Output:

```bash
      N E X T F L O W
      version 25.10.4 build 11173
      created 10-02-2026 15:17 UTC (16:17 CEST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

5. Verify where exactly Nextflow has been installed

```bash
which nextflow
```
Output:

```bash
/opt/anaconda3/envs/DNA/bin/nextflow
```

> [!IMPORTANT]
> Since Nextflow wasn't installed via conda, there's no Nextflow registration in conda's `DNA` metadata. Therefore, even though Nextflow won't appear in conda list of dependencies, it is physically there. Then, Nextflow is invisible to `conda list`, but don't worry, it was installed and functional!

> [!IMPORTANT]
> Now conda can work not only in `DNA` but also in `DNA2` Conda environments.


## 2. Installing Visual Studio Code (VS Code)

### Direct Download

1. Go to: <https://code.visualstudio.com/>

⚠️ IMPORTANT: If you have macOS version 11 (Big Sur) then go to this link:

<https://code.visualstudio.com/updates/v1_106>

to download the corresponding VS Code version. Choose '**Intel**' if your Mac was made before 2020.

2. Download macOS version

3. Drag to `/Applications`

4. Go to Launchpad or in Applications and click on the VS Code icon

Once VS Code is opened:

5. Install **Nextflow extension**

5.1. Go to the left panel and click on **Extensions**
5.2. Type in "Search Extensions in Marketplace" → "Nextflow"
5.3. Click on "Nextflow" (Nextflow language support)

⚠️ IMPORTANT: For those having macOS Big Sur, in VS Code, disable updates! Otherwise, VS Code will try to install/update the last version, provoking that the software can no longer work.

5.4 Open VS Code → Code → Preferences → Settings → Update: Mode → default

Change 'Mode' to: Update: Mode → none

## 3. Create folder structure

1. Go to ~/Genomics_cancer/data/

2. Create folder structure

mkdir -p SRR30536566_full_nf/{aligned,annotation,logs,qc,trimmed,variants}


## 4. Create nextflow script


1. See the folder structure of ~/Genomics_cancer

```code
Genomics_cancer/
├── data/
│   ├── SRR30536566/             # Used for step-by-step bash scripts (Part I & II)
│   │   ├── raw_fastq/           # `SRR30536566_1.fastq.gz`, `SRR30536566_2.fastq.gz`
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   ├── SRR30536566_full/             # Used by the unified bash script
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── logs/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   ├── SRR30536566_full_nf/   # Used by the unified nextflow script
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── logs/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   └── SRR30536566/   
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

2. Go to /scripts

3. Create the nextflow (.nf) file

touch 09_full_somatic_NGS_nextflow_script.nf

4. Open the .nf on VS Code

## 4. Check the amount of CPUs of your computer:

sysctl -n hw.ncpu

output:

4

## 5. Type the script




## 6. Run the .nf script from VS Code or Terminal

Go to: ~/Genomics_cancer/scripts

(DNA) $ nextflow run 09_full_somatic_NGS_nextflow_script.nf

 N E X T F L O W   ~  version 25.10.4

Launching `09_full_somatic_NGS_nextflow_script.nf` [voluminous_wilson] DSL2 - revision: 1d0af6be77

[-        ] FASTQC  -
[-        ] MULTIQC -
executor >  local (1)
[76/0f7269] FASTQC (SRR30536566) [  0%] 0 of 1
executor >  local (1)
[76/0f7269] FASTQC (SRR30536566) [  0%] 0 of 1
executor >  local (2)
[76/0f7269] FASTQC (SRR30536566) [100%] 1 of 1 ✔
executor >  local (2)
[76/0f7269] FASTQC (SRR30536566) [100%] 1 of 1 ✔
executor >  local (2)
[76/0f7269] FASTQC (SRR30536566) [100%] 1 of 1 ✔
[19/08e480] MULTIQC              [100%] 1 of 1



Update: Folder Structure

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
│       │     ├── Homo_sapiens_assembly38.fasta
│       │     ├── Homo_sapiens_assembly38.fasta.fai
│       │     ├── Homo_sapiens_assembly38.dict
│       │     ├── Homo_sapiens_assembly38.fasta.64.amb     
│       │     ├── Homo_sapiens_assembly38.fasta.64.ann     
│       │     ├── Homo_sapiens_assembly38.fasta.64.bwt     
│       │     ├── Homo_sapiens_assembly38.fasta.64.pac    
│       │     ├── Homo_sapiens_assembly38.fasta.64.sa     
│       │     └── Homo_sapiens_assembly38.fasta.64.alt    
│       │
│       ├── intervals/     
│       │     └── crc_panel_7genes_sorted.hg38.bed    
│       └── somatic_resources/    
│             ├── 1000g_pon.hg38.vcf.gz     
│             ├── af-only-gnomad.hg38.vcf.gz   
│             ├── 1000g_pon.hg38.vcf.gz.tbi    
│             └── af-only-gnomad.hg38.vcf.gz.tbi
│   
├── scripts/
│
└── logs/

```

(DNA2) Franos-MBP:Genomics_cancer Frano$ ls -lrth
total 8.0K
drwxr-xr-x  4 Frano staff  128 Jan  8 22:38 reference
drwxr-xr-x 14 Frano staff  448 Jan 15 23:28 logs
drwxr-xr-x  7 Frano staff  224 Mar  8 21:04 data
drwxr-xr-x 51 Frano staff 1.6K Mar 13 09:22 scripts


(DNA2) Franos-MBP:Genomics_cancer Frano$ ls -lrth data/
total 0
drwxr-xr-x 9 Frano staff 288 Jan  8 17:32 SRR30536566
drwxr-xr-x 9 Frano staff 288 Feb 27 23:43 SRR30536566_full_nf

(DNA2) Franos-MBP:Genomics_cancer Frano$ ls -lrth data/SRR30536566
total 0
drwxr-xr-x  5 Frano staff 160 Jan  8 18:11 raw_fastq

(DNA2) Franos-MBP:Genomics_cancer Frano$ ls -lrth data/SRR30536566_full_nf/
total 0
drwxr-xr-x 2 Frano staff 64 Feb 27 15:55 annotation
drwxr-xr-x 2 Frano staff 64 Feb 27 15:55 variants
drwxr-xr-x 3 Frano staff 96 Mar  5 21:47 aligned
drwxr-xr-x 3 Frano staff 96 Mar  6 00:01 trimmed
drwxr-xr-x 3 Frano staff 96 Mar  6 00:01 qc
drwxr-xr-x 3 Frano staff 96 Mar  6 00:01 logs

(DNA2) Franos-MBP:Genomics_cancer Frano$ 





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
│   │   │     ├── bwa_index.log  (optional)
│   │   │     └── SRR30536566_full_nf.flagstat.txt
│   │   ├── aligned/
│   │   │     ├── SRR30536566_full_nf.sorted.markdup.md.bam
│   │   │     ├── SRR30536566_full_nf.sorted.markdup.md.bam.bai
│   │   │     └── SRR30536566_full_nf.markdup.metrics.txt
│   │   ├── variants/
│   │   └── annotation/
│   │
├── reference/
│   └── GRCh38/
│       ├── fasta/
│       │     ├── Homo_sapiens_assembly38.fasta
│       │     ├── Homo_sapiens_assembly38.fasta.fai
│       │     ├── Homo_sapiens_assembly38.dict
│       │     ├── Homo_sapiens_assembly38.fasta.64.amb     
│       │     ├── Homo_sapiens_assembly38.fasta.64.ann     
│       │     ├── Homo_sapiens_assembly38.fasta.64.bwt     
│       │     ├── Homo_sapiens_assembly38.fasta.64.pac    
│       │     ├── Homo_sapiens_assembly38.fasta.64.sa     
│       │     └── Homo_sapiens_assembly38.fasta.64.alt    
│       │
│       ├── intervals/     
│       │     └── crc_panel_7genes_sorted.hg38.bed    
│       └── somatic_resources/    
│             ├── 1000g_pon.hg38.vcf.gz     
│             ├── af-only-gnomad.hg38.vcf.gz   
│             ├── 1000g_pon.hg38.vcf.gz.tbi    
│             └── af-only-gnomad.hg38.vcf.gz.tbi
│   
├── scripts/
│
└── logs/

```


## Conclusion

The bash script [09_full_somatic_NGS_bash_script.sh](bash_scripts/09_full_somatic_NGS_bash_script.sh) contains a pipeline that runs smoothly and outputting all expected files as the splitted bash scripts. Most importantly, this pipeline could also output the same three expected variants.

---

Go back to the top of  👉 [Part IV – Nextflow: Fully Automated Somatic DNA-NGS Pipeline](README_Part4_nextflow.md.md#part-iv--nextflow-fully-automated-somatic-dna-ngs-pipeline)

Visit the Bash script here 👉 [Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline](README_Part4_fullbash.md)

Go and see somatic NGS analysis in `DNA2` **samtools-updated** environment in 👉 [Part V: Pipeline maintenance and Environment Validation](README_Part5_DNA2_pipeline_update.md)

Jump to the first part of this tutorial 👉 [Part I – Preparation & setup](README_setup_Part1-3.md)

Go to the main page 👉 [Bash_pipeline_NGS](README.md)

