# Bash_pipeline_NGS

Welcome to my **DNA-NGS tutorial** 👋  

This repository shows you how to analyze **DNA sequencing datasets**  
using **Bash script pipelines**, starting from raw FASTQ data  
and running everything on a **home workstation or laptop**.

You can follow the tutorial step by step, or jump directly to the analysis part.

---

## Tutorial structure

### 1️⃣ Part I – Preparation & setup  
Learn how to prepare a clean and reproducible environment:
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
- Check potential artifacts and confirm annotated variants

➡️ **Go to analysis:**  
👉 [Part III – Variant Visualization](README_igv_Part3-3.md)

---

### 4️⃣ Part IV - Unified & Automated Pipeline (harmonized single-script somatic workflow)
- Bash script
- Nextflow

➡️ **Go to:**  
👉 [Part IV – Bash script: Fully Automated Somatic DNA-NGS Pipeline](README_fullbash.md)

## 🔮 Future extensions

This repository is designed to grow. Planned additions include:

- **Part V – Somatic analysis with matched panel of normals**
  - Additional datasets
  - Pipeline optimizations and best practices

- **Part VI – Germline analysis**
  - Additional datasets
  - Pipeline optimizations and best practices

---

## 🧬 Target audience

This tutorial is intended for:
- Bioinformatics students
- Life scientists learning NGS analysis
- Researchers who want a **transparent, Bash-only workflow**

---

## 📌 Notes

- The pipeline is optimized for **educational clarity**, not HPC clusters
- All steps are reproducible and runnable on a local machine
- Real public datasets from the **NCBI SRA** (e.g., SRR30536566) are used

---

Happy sequencing analysis!

---

## ⚠️ Disclaimer

This repository is intended for **educational and research purposes only**.  
It is not validated for clinical diagnostic use and should not be used for medical decision-making.

All analyses are performed on publicly available research datasets.

---

## 🔎 Third-Party Tools & Resources

This tutorial uses and displays output or screenshots generated from the following tools and databases:

- Integrative Genomics Viewer (IGV, Broad Institute)
- Genome Analysis Toolkit (GATK, Broad Institute)
- Ensembl Genome Browser
- OncoKB
- CIViC (Clinical Interpretation of Variants in Cancer)
- PanDrugs2
- FastQC
- MultiQC
- NCBI Sequence Read Archive (SRA)

All trademarks, software, and database contents belong to their respective owners.  
Screenshots and outputs are shown for educational and demonstration purposes only.

---

## 📜 License

© 2026 **bioinfo-frano**

This project is licensed under the **MIT License**. See the full license [here](LICENSE).
