# Bash_pipeline_NGS

Welcome to my **DNA-NGS tutorial** 👋  

This repository provides a step-by-step guide to analyzing **DNA sequencing datasets** using primarily **Bash script pipelines**, starting from raw FASTQ data and running the complete workflow on a **standard workstation or laptop**.

The tutorial can be followed sequentially for learning purposes or modularly depending on your goals.

This repository shows how to build and understand a complete somatic DNA-NGS analysis workflow, including:

- Data acquisition from public repositories
- Quality control and preprocessing
- Alignment and BAM processing
- Somatic variant calling (tumor-only)
- Filtering and biological post-processing
- Visualization and interpretation

The focus is on educational clarity, reproducibility, and transparency.
The workflow implements **tumor-only variant calling** without matched normal samples. However, one of the last part of this tutorial will be focused on **matched tumor-normal pair variant calling** analysis.

## 🔬 Workflow Overview
Overview of the tumor-only somatic DNA-NGS analysis pipeline implemented in this repository*.

<p align="center">
  <img src="images/Gemini_Generated_Image_2hszdp2hszdp2hsz.png" 
       alt="Tumor-only somatic DNA-NGS workflow" 
       width="65%">
</p>

- *Image generated in collaboration with Gemini (Google AI) via iterative prompting.*
---

## Tutorial structure

### 1️⃣ Part I – Preparation & setup  
Learn how to create a clean and reproducible environment:
- Folder structure
- Reference genome setup and integrity checks
- SRA data selection and download
- Conda environments and tool installation

➡️ **Start here:**  
👉 [Part I – Preparation & setup](README_setup_Part1-3.md)

---

### 2️⃣ Part II – Somatic analysis (Bash pipelines)  
Perform a **somatic DNA-NGS analysis** following GATK best practices:
- FASTQ processing and QC
- Alignment and BAM processing
- Somatic variant calling with **Mutect2**
- Variant filtering and annotation

➡️ **Go to analysis:**  
👉 [Part II – Somatic analysis](README_somatic_analysis_Part2-3.md)

---

### 3️⃣ Part III – Variant Visualization
- Learn how to visualize annotated variants in IGV
- Check potential artifacts and confirm annotated variant calls

➡️ **Go to analysis:**  
👉 [Part III – Variant Visualization](README_igv_Part3-3.md)

---

### 4️⃣ Part IV - Unified & Automated Pipeline (harmonized single-script somatic workflow)
A single, unified Bash script that runs the entire pipeline from end to end.

➡️ **Go to automated workflow:**  
👉 [Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline](README_Part4_fullbash.md)

---

### 5️⃣ Part V - Pipeline maintenance
Learn why it’s important to keep your pipeline up-to-date and ensure its reproducibility.

**Just click here**
👉 [Part V – Pipeline maintenance and Environment Validation](README_Part5_DNA2_pipeline_update.md)

---

### 6️⃣ Part VI - Nextflow Pipeline (harmonized single-script somatic workflow)
Understand what Nextflow is and how it compares to Bash pipelines. This section presents a fully automated Nextflow version of the somatic workflow.

➡️ **Go to automated workflow:**  
👉 [Part VI – Nextflow script: Fully Automated Somatic DNA-NGS Pipeline](README_Part6_Nextflow.md)

---

### 7️⃣ Part VII – Somatic analysis with matched tumor-normal pair
  - See the differences between tumor-only and matched tumor-normal somatic analysis
  - Get to know the advantages of this type of analysis compared to tumor-only

➡️ **Enter here:**  
👉 [Part VII – Matched Tumor-Normal Pair Somatic DNA-NGS Pipeline](README_Part7_tumor_normal.md)

---

## 🔮 Future extensions

This repository is designed to grow. Planned additions include:

- **Part VIII – Germline analysis**
  - New datasets and workflows for germline variant calling

---

## 🧬 Target audience

This tutorial is intended for:
- Bioinformatics students
- Life scientists learning NGS analysis
- Researchers who want a **transparent, Bash-only workflow** that can be run locally

---

## 📌 Notes

- The pipelines are optimized for **educational clarity** and designed to run on a standard workstation or laptop (not HPC clusters).
- All steps are reproducible and runnable on a local machine
- Real public datasets from the **NCBI SRA** (e.g., SRR30536566) are used

---

Happy sequencing analysis!

---

## 🧪 Tested Environment

- macOS (Intel)
- 8 GB RAM
- Free space: < 40GB
- macOS Big Sur 11.7.11 (Intel)
- Conda-based installation
- GATK 4.x

## ⚠️ Disclaimer

This repository is intended for **educational and research purposes only**.  
It is **not validated for clinical diagnostic use** and should not be used for medical decision-making.

All analyses are performed on publicly available research datasets.

---

## 🔎 Third-Party Tools & Resources

This tutorial uses and displays output or screenshots generated from the following tools and databases:

- Integrative Genomics Viewer (IGV, Broad Institute)
- Genome Analysis Toolkit (GATK, Broad Institute)
- Ensembl Genome Browser
- ClinVar
- OncoKB
- CIViC (Clinical Interpretation of Variants in Cancer)
- PanDrugs2
- FastQC
- MultiQC
- NCBI Sequence Read Archive (SRA)

All trademarks, software, and database contents belong to their respective owners.  
Screenshots and outputs are shown for educational and demonstration purposes only.

---

## 🤖 Acknowledgments

The Nextflow pipeline presented in **Part VI** was developed with the assistance of AI language models, specifically ChatGPT (OpenAI) and DeepSeek. These tools helped translate the original Bash workflow in **Part V** into a modular Nextflow script, and provided valuable guidance on debugging and optimisation. The Bash pipeline and all other parts of the tutorial were written by the author, with support from bioinformatics tool documentation, AI, and YouTube tutorials.

---

## 📜 License

© 2026 **bioinfo-frano**

This project is licensed under the **MIT License**. See the full license [here](LICENSE).
