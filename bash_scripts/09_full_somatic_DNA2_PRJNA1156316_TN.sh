#!/bin/bash

set -euo pipefail
set -o errtrace

trap 'echo "ERROR occurred at line $LINENO"; exit 1' ERR

LOG_PIPELINE="pipeline_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -i "$LOG_PIPELINE")
exec 2>&1

echo "======================================"
echo "Pipeline started: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"
echo "Conda env: $CONDA_DEFAULT_ENV"
echo "======================================"

echo "Using tools from:"
which samtools
which gatk
which bwa
which bcftools
which picard
which fastqc
which multiqc
which cutadapt

START_TIME=$(date +%s)

# ============================================================
# Configuration
# ============================================================
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

PROJECT_NAME="PRJNA1156316"      # <https://www.ebi.ac.uk/ena/browser/view/PRJNA1156316> <https://www.ncbi.nlm.nih.gov/sra/?term=SRR30536566> <https://www.ncbi.nlm.nih.gov/sra/?term=SRR30536541>
TUMOR_ID="SRR30536566"           # Sample Name: EIDR_55_tumor // Library Name: DMBEL-EIDR-071 // source_material_identifier: tumor   // tissue: colorectal cancer biopsy // sex: male
NORMAL_ID="SRR30536541"          # Sample Name: EIDR_55_blood // Library Name: DMBEL-EIDR-096 // source_material_identifier: blood 1 // tissue: blood                    // sex: male

SAMPLE_TUMOR="${TUMOR_ID}_tumor"       # Sample Name: EIDR_55_tumor
SAMPLE_NORMAL="${NORMAL_ID}_normal"    # Sample Name: EIDR_55_blood

PROJECT="$PROJECT_ROOT/data/$PROJECT_NAME"

THREADS=4

# NOTE:
# Variables that depend on SAMPLE must be defined inside the loop.
# Otherwise, Bash will fail with "unbound variable" when using 'set -u'.

# List of samples: Tumor first, Normal second
SAMPLES=("$SAMPLE_TUMOR" "$SAMPLE_NORMAL")

# Looping over each sample
for SAMPLE in "${SAMPLES[@]}"; do

  echo "======================================"
  echo "Processing sample for FastQC: $SAMPLE"
  echo "======================================"

  RAW_FASTQ_DIR="$PROJECT/$SAMPLE/raw_fastq"
  QC_DIR="$PROJECT/$SAMPLE/qc"

  mkdir -p "$QC_DIR/raw"

  echo " "
  echo "Starting FastQC for $SAMPLE..."
  echo " "

  fastqc \
    --threads "$THREADS" \
    --outdir "$QC_DIR/raw" \
    "$RAW_FASTQ_DIR"/*.fastq.gz

  echo " "
  echo "FastQC completed successfully for sample: $SAMPLE."
  echo " "

  echo "======================================"
  echo "Processing sample for MultiQC: $SAMPLE"
  echo "======================================"

  echo " "
  multiqc "$QC_DIR/raw" --outdir "$QC_DIR/raw"
  echo " "

  echo "MultiQC completed successfully for sample: $SAMPLE."

  echo " "


  echo "============================================================================"
  echo "Processing sample for Cutadapt trimming/filtering + full QC report: $SAMPLE"
  echo "============================================================================"

  echo " "

  RAW_FASTQ_DIR="$PROJECT/$SAMPLE/raw_fastq"
  QC_DIR="$PROJECT/$SAMPLE/qc"
  QC_DIR_TRIMMED="$PROJECT/$SAMPLE/qc/trimmed"

  #Cutadpat - trimming
  TRIM_DIR="$PROJECT/$SAMPLE/trimmed"
  LOG_DIR="$PROJECT/$SAMPLE/logs/"

  mkdir -p "$QC_DIR_TRIMMED"
  mkdir -p "$LOG_DIR"
  mkdir -p "$TRIM_DIR"

  cutadapt \
    -u 5 -u -5 \
    -U 5 -U -5 \
    -q 20,20 \
    -m 30 \
    -a A{10} \
    -A A{10} \
    -j "$THREADS" \
    --report=full \
    -o "$TRIM_DIR/${SAMPLE}_R1.trimmed.fastq.gz" \
    -p "$TRIM_DIR/${SAMPLE}_R2.trimmed.fastq.gz" \
    "$RAW_FASTQ_DIR/${SAMPLE}_1.fastq.gz" \
    "$RAW_FASTQ_DIR/${SAMPLE}_2.fastq.gz" \
    > "$LOG_DIR/cutadapt_${SAMPLE}.log"

  echo "Cutadapt completed successfully."


  # --- Run FastQC ---

  echo " "
  echo "FastQC starting for sample: $SAMPLE."
  echo " "

  fastqc \
    --threads "$THREADS" \
    --outdir "$QC_DIR_TRIMMED" \
    "$TRIM_DIR"/*.fastq.gz

  echo "Post trimming FastQC completed successfully."

  # --- Run MultiQC ---

  echo " "
  echo "MultiQC starting for sample: $SAMPLE."
  echo " "

  multiqc "$QC_DIR_TRIMMED" --outdir "$QC_DIR_TRIMMED"

  echo "Post trimming MultiQC completed successfully."

  echo " "

  echo "================================="
  echo "Step 0: Check / build BWA index"
  echo "================================="

  echo " "

  REF_DIR="$PROJECT_ROOT/reference/GRCh38/fasta"
  REF_FASTA="$REF_DIR/Homo_sapiens_assembly38.fasta"

  echo "Checking BWA index..."

  if [[ ! -f "${REF_FASTA}.64.bwt" && ! -f "${REF_FASTA}.bwt" ]]; then
    echo "BWA index not found for reference genome."
    echo "Building BWA index for reference genome..."
    bwa index "$REF_FASTA" 2> "$LOG_DIR/bwa_index.log"
    echo "BWA indexing completed."
  else
    echo "BWA index found. Skipping indexing of reference genome."
  fi

  echo " "

  echo "=========================================================================="
  echo "Step 1: Alignment: BWA-MEM with Read Groups → SAM. Sample: $SAMPLE"
  echo "=========================================================================="

  echo " "

  # Define alignment directory
  ALIGN_DIR="$PROJECT/$SAMPLE/aligned"
  mkdir -p "$ALIGN_DIR"

  # Assign correct SM (Sample name) tag depending on sample
  if [[ "$SAMPLE" == "$SAMPLE_TUMOR" ]]; then
    RG_SM="EIDR_55_tumor"
  elif [[ "$SAMPLE" == "$SAMPLE_NORMAL" ]]; then
    RG_SM="EIDR_55_blood"
  else
    echo "ERROR: Unknown sample $SAMPLE"
    exit 1
  fi

  # Read group information (REQUIRED by GATK)
  RG_ID="$SAMPLE"
  # RG_SM="EIDR_55_tumor || EIDR_55_blood"
  RG_LB="AMPLICON"
  RG_PL="ILLUMINA"
  RG_PU="HiSeq4000"

  echo "Running BWA-MEM alignment..."

  echo " "

  bwa mem \
    -t "$THREADS" \
    -R "@RG\tID:${RG_ID}\tSM:${RG_SM}\tLB:${RG_LB}\tPL:${RG_PL}\tPU:${RG_PU}" \
    "$REF_FASTA" \
    "$TRIM_DIR/${SAMPLE}_R1.trimmed.fastq.gz" \
    "$TRIM_DIR/${SAMPLE}_R2.trimmed.fastq.gz" \
    > "$ALIGN_DIR/${SAMPLE}.sam" \
    2> "$LOG_DIR/bwa_mem.log"

  echo "Alignment completed."

  echo " "

  echo "===================================================="
  echo " Step 2: Convert SAM → BAM: SAMPLE: $SAMPLE"
  echo "===================================================="

  echo " "

  if [[ ! -s "$ALIGN_DIR/${SAMPLE}.sam" ]]; then
    echo "ERROR: SAM file not created!" >&2
    exit 1
  fi

  echo " "

  echo "Converting SAM to BAM. Sample: $SAMPLE..."

  echo " "

  # Option: -S deprecated. Modern samtools auto-detects SAM/BAM
  samtools view \
    -@ "$THREADS" \
    -b \
    "$ALIGN_DIR/${SAMPLE}.sam" \
    > "$ALIGN_DIR/${SAMPLE}.bam"


  echo "================================================"
  echo "Step 3: Sort BAM (coordinate sort): SAMPLE: $SAMPLE"
  echo "================================================"

  echo " "

  if [[ ! -s "$ALIGN_DIR/${SAMPLE}.bam" ]]; then
    echo "ERROR: BAM file not created!" >&2
    exit 1
  fi

  echo "Sorting BAM. Sample: $SAMPLE..."

  echo " "

  samtools sort \
    -@ "$THREADS" \
    -o "$ALIGN_DIR/${SAMPLE}.sorted.bam" \
    "$ALIGN_DIR/${SAMPLE}.bam"

  # Check is sorted BAM file was created

  if [[ ! -s "$ALIGN_DIR/${SAMPLE}.sorted.bam" ]]; then
    echo "ERROR: Sorting failed!" >&2
    exit 1
  fi

  rm "$ALIGN_DIR/${SAMPLE}.sam" "$ALIGN_DIR/${SAMPLE}.bam"

  echo "Sorting completed."

  echo "============================================================================"
  echo "Step 4: Mark duplicates (AMPLICON-AWARE TAGGING OF DUPLICATES): SAMPLE: $SAMPLE"
  echo "============================================================================"

  echo " "

  if [[ ! -s "$ALIGN_DIR/${SAMPLE}.sorted.bam" ]]; then
    echo "ERROR: Sorted BAM not created!" >&2
    exit 1
  fi

  echo "Marking duplicates. Sample: $SAMPLE..."

  picard MarkDuplicates \
    INPUT="$ALIGN_DIR/${SAMPLE}.sorted.bam" \
    OUTPUT="$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam" \
    METRICS_FILE="$ALIGN_DIR/${SAMPLE}.markdup.metrics.txt" \
    CREATE_INDEX=false \
    REMOVE_DUPLICATES=false \
    TAG_DUPLICATE_SET_MEMBERS=true \
    VALIDATION_STRINGENCY=SILENT \
    2> "$LOG_DIR/markduplicates_${SAMPLE}.log"

  echo " "

  echo "Duplicate marking completed."

  echo " "

  echo "========================================================================================================"
  echo "Step 5: Add MD and NM tags (GATK robustness) + Index final BAM (REQUIRED for GATK): SAMPLE: $SAMPLE"
  echo "========================================================================================================"

  echo " "

  if [[ ! -s "$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam" ]]; then
    echo "ERROR: MarkDuplicates failed!" >&2
    exit 1
  fi

  # Define sorted markduplicates BAM files
  FINAL_BAM="$ALIGN_DIR/${SAMPLE}.sorted.markdup.md.bam"

  echo "Adding MD tags. Sample: $SAMPLE"

  samtools calmd \
    -b \
    "$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam" \
    "$REF_FASTA" \
    > "$FINAL_BAM"

  if [[ ! -s "$FINAL_BAM" ]]; then
    echo "ERROR: Final BAM not created by samtools calmd!" >&2
    exit 1
  fi

  echo "Indexing final BAM. Sample: $SAMPLE"

  samtools index "$FINAL_BAM"

  echo "BAM indexing completed."

  echo " "

  echo "=========================================="
  echo "Step 6: Alignment statistics: $SAMPLE"
  echo "=========================================="

  echo " "

  echo "Generating alignment statistics. Sample: $SAMPLE"

  samtools flagstat \
    "$FINAL_BAM" \
    > "$LOG_DIR/${SAMPLE}.flagstat.txt"

  samtools idxstats \
    "$FINAL_BAM" \
    > "$LOG_DIR/${SAMPLE}.idxstats.txt"

  echo " "

  echo "========================================================================"
  echo "Step 7: MultiQC (duplicates + alignment metrics). SAMPLE: $SAMPLE"
  echo "========================================================================"

  echo " "

  # Create subfolder in ~/SAMPLE/qc

  mkdir -p "$QC_DIR/md_flagstat"

  echo "Running MultiQC. Sample: $SAMPLE"

  multiqc \
    "$ALIGN_DIR" \
    "$LOG_DIR" \
    --outdir "$QC_DIR/md_flagstat"

  echo "MultiQC of MarkDuplicates & flagstat & Cutadapt done!"

  echo " "

  echo "=================================================================="
  echo "Step 8: Cleanup intermediate files (Optional). SAMPLE: $SAMPLE"
  echo "=================================================================="

  echo "Cleaning up intermediate files..."

  rm -f \
    "$ALIGN_DIR/${SAMPLE}.sorted.bam" \
    "$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam"

  echo "Cleanup completed."

  echo " "
  echo "Alignment and BAM preprocessing completed successfully."
  echo " "

done


# =========================================
# Somatic variant calling: Configuration
# =========================================

echo " "

JAVA_MEM="-Xmx6g"

VARIANT_DIR="$PROJECT/variants"
LOG_DIR_PROJECT="$PROJECT/logs"

REF_DIR="$PROJECT_ROOT/reference/GRCh38/fasta"
REF_FASTA="$REF_DIR/Homo_sapiens_assembly38.fasta"
REF_FASTA_DICT="$REF_DIR/Homo_sapiens_assembly38.dict"
INTERVALS="$PROJECT_ROOT/reference/GRCh38/intervals/crc_panel_7genes_sorted.hg38.bed"
SOMATIC_RESOURCES="$PROJECT_ROOT/reference/GRCh38/somatic_resources"
PON="$SOMATIC_RESOURCES/1000g_pon.hg38.vcf.gz"
GNOMAD="$SOMATIC_RESOURCES/af-only-gnomad.hg38.vcf.gz"

ALIGN_DIR_TUMOR="$PROJECT/$SAMPLE_TUMOR/aligned"
ALIGN_DIR_NORMAL="$PROJECT/$SAMPLE_NORMAL/aligned"

TUMOR_MD_BAM="$ALIGN_DIR_TUMOR/${SAMPLE_TUMOR}.sorted.markdup.md.bam"
NORMAL_MD_BAM="$ALIGN_DIR_NORMAL/${SAMPLE_NORMAL}.sorted.markdup.md.bam"

# Define alignment directory
mkdir -p "$VARIANT_DIR"
mkdir -p "$LOG_DIR_PROJECT"

# ============================================
# Sanity checks (fail early, fail clearly)
# ============================================

echo "Checking required inputs for Mutect2..."

for file in \
  "$REF_FASTA" \
  "$REF_FASTA_DICT" \
  "$REF_FASTA.fai" \
  "$TUMOR_MD_BAM" \
  "$TUMOR_MD_BAM.bai" \
  "$NORMAL_MD_BAM" \
  "$NORMAL_MD_BAM.bai" \
  "$INTERVALS" \
  "$PON" \
  "$PON.tbi" \
  "$GNOMAD" \
  "$GNOMAD.tbi"
do
  [[ -s "$file" ]] || { echo "ERROR: Missing required file: $file"; exit 1; }
done

echo " "

echo "All required input files found."

echo " "

echo "======================================================================================"
echo  "Running Mutect2 (matched tumor-normal panel of 7 genes). Project: $PROJECT_NAME"
echo "======================================================================================"

TUMOR_SM="EIDR_55_tumor"
NORMAL_SM="EIDR_55_blood"

gatk --java-options "$JAVA_MEM" Mutect2 \
  -R "$REF_FASTA" \
  -I "$TUMOR_MD_BAM" \
  -I "$NORMAL_MD_BAM" \
  --tumor-sample "$TUMOR_SM" \
  --normal-sample "$NORMAL_SM" \
  --panel-of-normals "$PON" \
  --germline-resource "$GNOMAD" \
  -L "$INTERVALS" \
  --f1r2-tar-gz "$VARIANT_DIR/${PROJECT_NAME}.f1r2.tar.gz" \
  -O "$VARIANT_DIR/${PROJECT_NAME}_tumor_normal.unfiltered.vcf.gz" \
  > "$LOG_DIR_PROJECT/mutect2.stdout.log" \
  2> "$LOG_DIR_PROJECT/mutect2.stderr.log"

echo " "
echo "Mutect2 completed successfully."
echo "Unfiltered VCF written to: $VARIANT_DIR/${PROJECT_NAME}_tumor_normal.unfiltered.vcf.gz"

## --af-of-alleles-not-in-resource 0.0000025 \ --> A must-use in tumor-only, better remove it from tumor-normal analysis


echo "========================================="
echo "LearnReadOrientationModel: Configuration"
echo "========================================="

VARIANT_DIR="$PROJECT/variants"
LOG_DIR_PROJECT="$PROJECT/logs"

F1R2_TAR="$VARIANT_DIR/${PROJECT_NAME}.f1r2.tar.gz"
ORIENTATION_MODEL="$VARIANT_DIR/${PROJECT_NAME}.read-orientation-model.tar.gz"

echo " "
echo "Checking required inputs for LearnReadOrientationModel..."
echo " "
echo "Running sanity checks..."

# ===============================================================
# Sanity checks: LearnReadOrientationModel // GetPileupSummaries
# ===============================================================

# Sanity checks
for file in \
  "$F1R2_TAR" \
  "$TUMOR_MD_BAM" \
  "$TUMOR_MD_BAM.bai" \
  "$NORMAL_MD_BAM" \
  "$NORMAL_MD_BAM.bai" \
  "$REF_FASTA" \
  "$REF_FASTA.fai" \
  "$REF_FASTA_DICT" \
  "$GNOMAD" \
  "$GNOMAD.tbi"
do
  [[ -s "$file" ]] || { echo "ERROR: Missing required file: $file"; exit 1; }
done
echo " "
echo "All required files found."

echo "============================================================"
echo  "LearnReadOrientationModel. Project: $PROJECT_NAME"
echo "============================================================"

echo " "

gatk --java-options "$JAVA_MEM" LearnReadOrientationModel \
  -I "$F1R2_TAR" \
  -O "$ORIENTATION_MODEL" \
  2> "$LOG_DIR_PROJECT/${PROJECT_NAME}_learn_read_orientation_model.log"

# If $ORIENTATION_MODEL output is empty -> "Error"
[[ -s "$ORIENTATION_MODEL" ]] || { echo "Orientation model failed"; exit 1; }

echo "LearnReadOrientationModel completed successfully."
echo "Orientation model written to: $ORIENTATION_MODEL"

echo " "

echo "============================================================"
echo "GetPileupSummaries. Project: $PROJECT_NAME"
echo "============================================================"

echo " "

PILEUP_TUMOR="$VARIANT_DIR/${PROJECT_NAME}_tumor.pileups.table"
PILEUP_NORMAL="$VARIANT_DIR/${PROJECT_NAME}_normal.pileups.table"

gatk --java-options "$JAVA_MEM" GetPileupSummaries \
  -R "$REF_FASTA" \
  -I "$TUMOR_MD_BAM" \
  -V "$GNOMAD" \
  -L "$INTERVALS" \
  -O "$PILEUP_TUMOR" \
  2> "$LOG_DIR_PROJECT/${PROJECT_NAME}_${SAMPLE_TUMOR}_get_pileup_summaries.log"

echo "GetPileupSummaries completed for: $PILEUP_TUMOR."

# To get to know the amount of informative SNPs. If number of sites < 10–20, contamination estimate is statistically weak.

echo " "

gatk --java-options "$JAVA_MEM" GetPileupSummaries \
  -R "$REF_FASTA" \
  -I "$NORMAL_MD_BAM" \
  -V "$GNOMAD" \
  -L "$INTERVALS" \
  -O "$PILEUP_NORMAL" \
  2> "$LOG_DIR_PROJECT/${PROJECT_NAME}_${SAMPLE_NORMAL}_get_pileup_summaries.log"

echo "GetPileupSummaries completed for: $PILEUP_NORMAL."

echo " "

SITE_COUNT_TUMOR=$(grep -v "^#" "$PILEUP_TUMOR" | wc -l)
SITE_COUNT_NORMAL=$(grep -v "^#" "$PILEUP_NORMAL" | wc -l)

echo "Tumor SNP sites: $SITE_COUNT_TUMOR"
echo "Normal SNP sites: $SITE_COUNT_NORMAL"

if [[ "$SITE_COUNT_TUMOR" -lt 10 ]]; then
  echo "WARNING: Very few informative SNPs ($SITE_COUNT_TUMOR). Contamination estimate may be unreliable."
fi

echo "CalculateContamination will start after sanity checks"
echo " "
# ============================================================
# Sanity checks
# ============================================================

echo "Running sanity checks..."

for file in "$PILEUP_TUMOR" "$PILEUP_NORMAL"
do
  [[ -s "$file" ]] || { echo "ERROR: Missing $file"; exit 1; }
done

echo "All required files found."

echo " "
echo "============================================================"
echo "CalculateContamination. Project: $PROJECT_NAME"
echo "============================================================"

CONTAM_TABLE="$VARIANT_DIR/${PROJECT_NAME}.contamination.table"
SEGMENTS_TABLE="$VARIANT_DIR/${PROJECT_NAME}.segments.table"

echo "Starting CalculateContamination..."

gatk --java-options "$JAVA_MEM" CalculateContamination \
  -I "$PILEUP_TUMOR" \
  -matched "$PILEUP_NORMAL" \
  -O "$CONTAM_TABLE" \
  --tumor-segmentation "$SEGMENTS_TABLE" \
  2> "$LOG_DIR_PROJECT/${PROJECT_NAME}_calculate_contamination.log"

echo "CalculateContamination completed."
echo "Contamination table written to: $CONTAM_TABLE"

# If $CONTAM_TABLE output is empty -> "Error"
[[ -s "$CONTAM_TABLE" ]] || { echo "Contamination estimation failed"; exit 1; }


# ============================================================
# Sanity checks: FilterMutectCalls
# ============================================================

echo " "
echo "Variant filtering will start after sanity checks"

REF_DIR="$PROJECT_ROOT/reference/GRCh38/fasta"
REF_FASTA="$REF_DIR/Homo_sapiens_assembly38.fasta"
LOG_DIR_PROJECT="$PROJECT/logs"
ORIENTATION_MODEL="$VARIANT_DIR/${PROJECT_NAME}.read-orientation-model.tar.gz"
CONTAM_TABLE="$VARIANT_DIR/${PROJECT_NAME}.contamination.table"
SEGMENTS_TABLE="$VARIANT_DIR/${PROJECT_NAME}.segments.table"
UNFILTERED_VCF="$VARIANT_DIR/${PROJECT_NAME}_tumor_normal.unfiltered.vcf.gz"
FILTERED_VCF="$VARIANT_DIR/${PROJECT_NAME}.filtered.vcf.gz"

echo " "
echo "Running sanity checks..."

for file in \
  "$REF_FASTA" \
  "$UNFILTERED_VCF" \
  "$UNFILTERED_VCF.tbi" \
  "$ORIENTATION_MODEL" \
  "$CONTAM_TABLE"
do
  [[ -s "$file" ]] || { echo "ERROR: Missing required file: $file"; exit 1; }
done

echo " "
echo "All required files found."
echo " "

echo "============================================================"
echo "FilterMutectCalls. Project: $PROJECT_NAME"
echo "============================================================"

echo " "

gatk --java-options "$JAVA_MEM" FilterMutectCalls \
  -R "$REF_FASTA" \
  -V "$UNFILTERED_VCF" \
  --contamination-table "$CONTAM_TABLE" \
  --orientation-bias-artifact-priors "$ORIENTATION_MODEL" \
  --tumor-segmentation "$SEGMENTS_TABLE" \
  -O "$FILTERED_VCF" \
  2> "$LOG_DIR_PROJECT/${PROJECT_NAME}_filter_mutect_calls.log"

# Check if $FILTERED_VCF output exists and is not empty -> otherwise "Error"
[[ -s "$FILTERED_VCF" ]] || { echo "FilterMutectCalls failed"; exit 1; }

echo "FilterMutectCalls completed successfully."
echo "Final filtered VCF: $FILTERED_VCF"

echo " "

echo "============================================================"
echo "Post-Filtering. Project: $PROJECT_NAME"
echo "============================================================"

# ============================================================
# Extract PASS variants only
# ============================================================

TUMOR_SM="EIDR_55_tumor"
NORMAL_SM="EIDR_55_blood"
FILTERED_VCF="$VARIANT_DIR/${PROJECT_NAME}.filtered.vcf.gz"
PASS_VCF="$VARIANT_DIR/${PROJECT_NAME}.filtered.PASS.vcf.gz"
LOG_FILE="$LOG_DIR_PROJECT/${PROJECT_NAME}.postfilter.log"
POSTFILTER_VCF="$VARIANT_DIR/${PROJECT_NAME}.postfiltered.vcf.gz"
SUMMARY_TXT="$VARIANT_DIR/${PROJECT_NAME}.postfilter_summary.txt"

echo "Extracting PASS variants..."

bcftools view -f PASS "$FILTERED_VCF" -Oz -o "$PASS_VCF"       # "-f PASS" only shows variants where the FILTER column == PASS. It will return only high-confidence calls.
bcftools index -t "$PASS_VCF"                                  # "-t" when indexing as ".tbi", otherwise ".csi" by default.

# Count PASS variants
PASS_COUNT=$(bcftools view -H "$PASS_VCF" | wc -l)             # "-H" hides header. It reads a file that already contains only PASS variants (from the '-f PASS code')
echo "Number of PASS variants: $PASS_COUNT"

echo "PASS-only VCF written to: $PASS_VCF"

echo " "
bcftools query -l $FILTERED_VCF

# Output:
# EIDR_55_blood   → index 0   → NORMAL = sample 0
# EIDR_55_tumor   → index 1   → TUMOR = sample 1
# The order is assigned by GATK alphabetically by sample name

echo " "
echo "Post variant filtering will start after sanity checks"
echo " "

# ============================================================
# Thresholds (amplicon tumor-normal)
# ============================================================

MIN_DP=200          # total depth
MIN_AD_ALT=10       # ALT read count
MIN_VAF=0.02        # 2%

# ============================================================
# Logging
# ============================================================

exec > >(tee -a "$LOG_FILE") 2>&1

echo "Starting post-filtering. Sample: $PROJECT_NAME"

# ============================================================
# Sanity checks
# ============================================================

for file in "$FILTERED_VCF" "$FILTERED_VCF.tbi"
do
  [[ -f "$file" ]] || { echo "ERROR: Missing $file"; exit 1; }
done


# ============================================================
# Apply hard filters
# ============================================================
echo " "
echo "Applying post-filter thresholds:"
echo "  DP >= $MIN_DP"
echo "  ALT reads (AD[1]) >= $MIN_AD_ALT"
echo "  VAF >= $MIN_VAF"

# bcftools view -s "$TUMOR_SM" -f PASS -Ou "$FILTERED_VCF" -> Extracting only the tumor sample from $FILTERED_VCF, so the resulting $POSTFILTER_VCF has only one sample: EIDR_55_tumor
# Keep both tumor and normal, tumor first
bcftools view -s "$TUMOR_SM","$NORMAL_SM" -f PASS -Ou "$FILTERED_VCF" | \
bcftools filter \
  -i "FORMAT/DP >= ${MIN_DP} && FORMAT/AD[0:1] >= ${MIN_AD_ALT} && FORMAT/AF >= ${MIN_VAF}" \
  -Oz -o "$POSTFILTER_VCF"

# Only tumor sample is present
# FORMAT/AD[0:1] correctly refers to sample 0 (tumor), ALT allele
echo " "
bcftools query -l "$FILTERED_VCF"
echo " "
bcftools query -l "$POSTFILTER_VCF"
# Should return:
# EIDR_55_tumor
# EIDR_55_blood

# ============================================================
# Index the post-filtered VCF (required for IGV) → .tbi
# ============================================================
echo " "
echo "Indexing post-filtered VCF"

bcftools index -t "$POSTFILTER_VCF"                           # Option '-t' → .tbi. Without any option → .csi (default). Both are valid index.

# Sanity check: ensure index was created
if [[ ! -f "${POSTFILTER_VCF}.tbi" ]]; then
  echo "ERROR: Tabix index (.tbi) was not created"
  exit 1
fi

# ============================================================
# Variant counts
# ============================================================
echo " "

N_VARIANTS=$(bcftools view -H "$POSTFILTER_VCF" | wc -l)

if [[ "$N_VARIANTS" -eq 0 ]]; then
  echo "WARNING: 0 variants passed post-filtering"
else
  echo "Variants retained after post-filtering: $N_VARIANTS"
fi

# ============================================================
# Summary file
# ============================================================

{
  echo "Post-filter summary"
  echo "========================"
  echo "Project: $PROJECT_NAME"
  echo "Date: $(date)"
  echo ""
  echo "Library type: Amplicon (PCR)"
  echo "Sequencing: Tumor-Normal (paired)"
  echo "Post-filtering applied on tumor sample only, normal sample retained for paired analysis"
  echo ""
  echo "Tumor sample: $TUMOR_SM"
  echo "Normal sample: $NORMAL_SM"
  echo ""
  echo "Input VCF: $FILTERED_VCF"
  echo ""
  echo "Thresholds:"
  echo "  DP >= $MIN_DP"
  echo "  ALT reads >= $MIN_AD_ALT"
  echo "  VAF >= $MIN_VAF"
  echo ""
  echo "PASS variants before post-filtering: $PASS_COUNT"
  echo "Variants retained: $N_VARIANTS"
  echo ""
  PCT=$(awk -v a="$N_VARIANTS" -v b="$PASS_COUNT" 'BEGIN {if (b>0) printf "%.2f", (a/b)*100; else print "0"}')
  echo "Retention rate: ${PCT}%"

} > "$SUMMARY_TXT"

echo " "
echo "[$(date)] $SAMPLE_TUMOR/$TUMOR_SM post-filtering completed successfully"
echo " "

# ========================
# TOTAL PIPELINE RUNTIME
# ========================

END_TIME_TOTAL=$(date +%s)
TOTAL_ELAPSED=$((END_TIME_TOTAL - START_TIME))

echo "=========================================="
echo "FULL PIPELINE COMPLETED SUCCESSFULLY"
echo "Total runtime: ${TOTAL_ELAPSED} seconds"
echo "Total runtime: $((TOTAL_ELAPSED / 60)) minutes"
echo "Total runtime: $((TOTAL_ELAPSED / 3600)) hours"
echo "=========================================="
