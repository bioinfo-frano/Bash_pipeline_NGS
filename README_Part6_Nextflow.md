# Part VI – Nextflow Pipeline: Fully Automated Somatic DNA-NGS Pipeline (Single Nextflow Script)

In [Part IV – Bash script](README_Part4_fullbash.md) we introduced a single bash script that performs the complete NGS somatic analysis of the 7-gene amplicon dataset `SRR30536566`. This time we present a **Nextflow pipeline**  that accomplishes the same analysis.

## What is Nextflow?

**Nextflow** is an open-source, data-driven workflow management system designed for building scalable, portable, and reproducible computational pipelines, especially in bioinformatics. The scripting language (based on Groovy) has a learning curve, and its concepts (channels, processes, operators) may feel unfamiliar at first, but once mastered, it becomes extremely powerful.
Nextflow can run workflows on your local machine, HPC clusters, or cloud platforms with minimal changes.

## Advantages and Disadvantages: Nextflow vs. Bash

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

- [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)

- [Maxwell Cluster – Nextflow guide](https://wiki.desy.de/maxwell/documentation/workflows/nextflow/)

- [nf-core community](https://nf-co.re/) – "*A global community collaborating to build open-source Nextflow components and pipelines*" 

---

## Nextflow analysis of dataset `SRR30536566`

Nextflow is a **Java-based workflow manager** that is updated frequently. The official installation method (`curl -s https://get.nextflow.io | bash`) is recommended for several reasons:

1. **Latest version immediately** – The curl method always downloads the most recent stable release directly from the Nextflow website. Conda packages may lag behind by weeks or months, meaning you could miss important features, performance improvements, or bug fixes.

2. **Simplicity and control** – Nextflow is a single executable file (a Bash wrapper that launches Java). Installing it manually gives you full control over which version you use and where it resides. Updating is as simple as re-running the curl command or replacing the executable.

3. **No dependency on Conda's Python ecosystem** – Nextflow only requires Java (version 11 or later, with 17+ recommended). It does not rely on Python, R, or any Conda-managed libraries. Installing it via Conda adds unnecessary metadata and couples it to a specific environment, which can complicate maintenance.

4. **Official recommendation** – The Nextflow developers explicitly recommend the manual installation method for production use (see the [official documentation](https://www.nextflow.io/docs/latest/install.html#installation)). Conda packages are community-maintained and may not always be up-to-date.

5. **Portability** – Because Nextflow is a single file, you can easily move it between Conda environments (as we will do) or even keep a shared copy in a system-wide location (like `/usr/local/bin`). This flexibility is lost when using a Conda package.

Therefore, Nextflow will be installed manually (not via Conda).
The executable will be placed inside the `DNA` environment's `bin` directory so that it is available when the environment is active.

---

## 1.  Install Nextflow (manually) in the `DNA` environment

You can use the same mechanism to install it in `DNA2` later.

**In Terminal**

1. Activate `DNA` environment and check that Java version >= 17 is available. 

```bash
conda activate DNA
java -version
```
Output:

```bash
openjdk version "17.0.17" 2025-10-21 LTS
```

2. Download and install Nextflow using the official script:

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

>[!IMPORTANT] This command is not installing Nextflow via conda, just downloading an executable file called `nextflow` in your current directory. That means that by typing "conda list", you won't find "nextflow".
 

3. Now move it into your Conda environment’s `bin` folder and make it executable:

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

5. Confirm location of Nextflow:

```bash
which nextflow
```
Output:

```bash
/opt/anaconda3/envs/DNA/bin/nextflow
```

> [!IMPORTANT]
> Because Nextflow was installed manually (via `curl`), it does not appear in conda list. However, it is physically present in the environment's bin directory and fully functional, so don't worry, it was installed and functional!

> [!IMPORTANT]
> Nextflow is now installed in the `DNA` environment, but it can also be used with other environments (like `DNA2`) because Nextflow itself is environment‑agnostic. The actual tools used by the pipeline are specified via the `conda "DNA2"` directive inside the script. You can either keep Nextflow in `DNA` or move it to `DNA2` – both work fine.

## 2. (Optional) Move Nextflow to the `DNA2` environment

If you prefer to keep Nextflow inside the `DNA2` environment (which contains the updated samtools), you can easily move it:

```bash
conda activate DNA2
mv /opt/anaconda3/envs/DNA/bin/nextflow $CONDA_PREFIX/bin/
chmod +x $CONDA_PREFIX/bin/nextflow
which nextflow   # should point to DNA2/bin/nextflow
```

## 3. (Optional) Setting up Visual Studio Code for Nextflow script editing

### Direct Download

1. Go to: <https://code.visualstudio.com/>

>[!IMPORTANT]  If you have macOS version 11 (Big Sur) then go to this link:
>
><https://code.visualstudio.com/updates/v1_106>
>
>to download the corresponding VS Code version. Choose '**Intel**' if your Mac was made before 2020.

2. Download macOS version

3. Drag to `/Applications`

4. Go to Launchpad or in Applications and click on the VS Code icon

Once VS Code is opened:

5. Install **Nextflow extension**

5.1. Go to the left panel and click on **Extensions**
5.2. Type in "Search Extensions in Marketplace" → "Nextflow"
5.3. Click on "Nextflow" (Nextflow language support)

>[!IMPORTANT] For those having macOS Big Sur, it's recommended to **disable updates** in VS Code. Otherwise, VS Code will try to install/update the last version, otherwise VS Code may attempt to install a newer, incompatible version.

6. Disabling updates:

I. Open VS Code
II. Go to **Code** → **Preferences** → **Settings**
III. Search for "update mode"
IV. Change **Update: Mode** from `default` to `none`


## 4. Create folder structure

1. Go to ~/Genomics_cancer/data/

2. Create folder structure

mkdir -p SRR30536566_full_nf/{aligned,annotation,logs,qc,trimmed,variants}


## 5. Create nextflow script


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
│   ├── SRR30536566_full/        # Used by the unified bash script
│   │   ├── qc/
│   │   ├── trimmed/
│   │   ├── logs/
│   │   ├── aligned/
│   │   ├── variants/
│   │   └── annotation/
│   │
│   ├── SRR30536566_full_nf/    # Used by the unified nextflow script
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

```bash
touch 09_full_somatic_Nextflow.nf
```
4. Open the .nf on VS Code

## 6. Check the amount of CPUs of your computer:

```bash
sysctl -n hw.ncpu
```
output: `4`

## 7. The full Nextflow script

👉 [09_full_somatic_home_nextflow.nf](nextflow_scripts/09_full_somatic_home_nextflow.nf)

1. Download the script and place it in `~/Genomics_cancer/scripts`

2. Open it and manually replace  `home_dir = System.getenv('HOME')` for your own path to `/Genomics_cancer/scripts`.

For example: `home_dir = /path/to/your`

Then the script will look for files in `/path/to/your/Genomics_cancer/...`

>[!IMPORTANT] Make sure the folder structure exactly matches the paths used.

## 8. Run the `.nf` script from VS Code or Terminal

1. Activate `DNA2` and go to `~/Genomics_cancer/scripts`

2. Run:

```bash
(DNA2) scripts $ nextflow run 09_full_somatic_home_nextflow.nf

 N E X T F L O W   ~  version 25.10.4

Launching `script_test8_full.nf` [ridiculous_rubens] DSL2 - revision: 6410c369f4

executor >  local (13)
[24/603d41] process > FASTQC (SRR30536566)         [100%] 1 of 1 ✔
[43/301f80] process > MULTIQC_RAW                  [100%] 1 of 1 ✔
[38/44130b] process > CUTADAPT (SRR30536566)       [100%] 1 of 1 ✔
[ef/d597a2] process > FASTQC_TRIMMED (SRR30536566) [100%] 1 of 1 ✔
[a8/62b73e] process > MULTIQC_TRIMMED              [100%] 1 of 1 ✔
[d9/c9d855] process > ALIGNMENT (SRR30536566)      [100%] 1 of 1 ✔
[7e/15bc70] process > MULTIQC_ALIGNMENT (1)        [100%] 1 of 1 ✔
[a5/0f473b] process > MUTECT2 (1)                  [100%] 1 of 1 ✔
[3c/121e60] process > LEARN_READ_ORIENTATION (1)   [100%] 1 of 1 ✔
[f4/25df68] process > GET_PILEUP_SUMMARIES (1)     [100%] 1 of 1 ✔
[f2/8e6f2b] process > CALCULATE_CONTAMINATION (1)  [100%] 1 of 1 ✔
[de/601d26] process > FILTER_MUTECT_CALLS (1)      [100%] 1 of 1 ✔
[0e/340e22] process > POSTFILTER_VARIANTS (1)      [100%] 1 of 1 ✔
Completed at: 14-Mar-2026 10:59:53
Duration    : 28m 10s
CPU hours   : 1.8
Succeeded   : 13
```

## 9. Outputted files (qc, trimmed, logs, aligned, variants) Folder Structure

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
│   └── 09_full_somatic_home_Nextflow.nf 
│
└── logs/


```


## Conclusion

The nextflow script [09_full_somatic_home_nextflow.nf](nextflow_scripts/09_full_somatic_home_nextflow.nf) contains a single pipeline that runs smoothly end-to-end the somatic analysis, producing all expected files as the individual Bash scripts in [Part II – Somatic analysis](README_somatic_analysis_Part2-3.md#part-ii--somatic-analysis-bash-pipelines). The Nextflow pipeline could also output the same three expected variants.

---

Go back to the top of  👉 [Part VI – Nextflow: Fully Automated Somatic DNA-NGS Pipeline](README_Part6_nextflow.md.md#part-vi--nextflow-fully-automated-somatic-dna-ngs-pipeline)

Visit the Bash script here 👉 [Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline](README_Part4_fullbash.md)

Go and see somatic NGS analysis in `DNA2` **samtools-updated** environment in 👉 [Part V: Pipeline maintenance and Environment Validation](README_Part5_DNA2_pipeline_update.md)

Jump to the first part of this tutorial 👉 [Part I – Preparation & setup](README_setup_Part1-3.md)

Go to the main page 👉 [Bash_pipeline_NGS](README.md)

