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

echo "Starting FastQC..."

# ============================================================
# Configuration
# ============================================================
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="SRR30536566"
SAMPLE_FULL="${SAMPLE}_full_DNA2"
TRIM_DIR="$PROJECT_ROOT/data/$SAMPLE_FULL/trimmed"

RAW_FASTQ_DIR="$PROJECT_ROOT/data/$SAMPLE/raw_fastq"
QC_DIR="$PROJECT_ROOT/data/$SAMPLE_FULL/qc"
LOG_DIR="$PROJECT_ROOT/data/$SAMPLE_FULL/logs/"

THREADS=4

echo "PROJECT_ROOT: $PROJECT_ROOT"
echo "SAMPLE: $SAMPLE"
echo "SAMPLE_FULL: $SAMPLE_FULL"

REF_DIR="$PROJECT_ROOT/reference/GRCh38/fasta"
REF_FASTA="$REF_DIR/Homo_sapiens_assembly38.fasta"
REF_FASTA_DICT="$REF_DIR/Homo_sapiens_assembly38.dict"
INTERVALS="$PROJECT_ROOT/reference/GRCh38/intervals/crc_panel_7genes_sorted.hg38.bed"

ALIGN_DIR="$PROJECT_ROOT/data/$SAMPLE_FULL/aligned"

FINAL_BAM="$ALIGN_DIR/${SAMPLE_FULL}.sorted.markdup.md.bam"

# Read group information (REQUIRED by GATK)
RG_ID="SRR30536566"
RG_SM="DMBEL-EIDR-071"
RG_LB="AMPLICON"
RG_PL="ILLUMINA"
RG_PU="HiSeq4000"


SOMATIC_RESOURCES="$PROJECT_ROOT/reference/GRCh38/somatic_resources"
PON="$SOMATIC_RESOURCES/1000g_pon.hg38.vcf.gz"
GNOMAD="$SOMATIC_RESOURCES/af-only-gnomad.hg38.vcf.gz"

BAM_DIR="$PROJECT_ROOT/data/$SAMPLE_FULL/aligned"
INPUT_MD_BAM="$BAM_DIR/${SAMPLE_FULL}.sorted.markdup.md.bam"
TUMOR_SM="DMBEL-EIDR-071"                                                           # biological (tumor) sample name

JAVA_MEM="-Xmx6g"

VARIANT_DIR="$PROJECT_ROOT/data/$SAMPLE_FULL/variants"
F1R2_TAR="$VARIANT_DIR/${SAMPLE_FULL}.f1r2.tar.gz"
ORIENTATION_MODEL="$VARIANT_DIR/${SAMPLE_FULL}.read-orientation-model.tar.gz"

PILEUP_TABLE="$VARIANT_DIR/${SAMPLE_FULL}.pileups.table"
CONTAM_TABLE="$VARIANT_DIR/${SAMPLE_FULL}.contamination.table"

UNFILTERED_VCF="$VARIANT_DIR/${SAMPLE_FULL}.unfiltered.vcf.gz"
FILTERED_VCF="$VARIANT_DIR/${SAMPLE_FULL}.filtered.vcf.gz"
PASS_VCF="$VARIANT_DIR/${SAMPLE_FULL}.filtered.PASS.vcf.gz"

POSTFILTER_VCF="$VARIANT_DIR/${SAMPLE_FULL}.postfiltered.vcf.gz"
SUMMARY_TXT="$VARIANT_DIR/${SAMPLE_FULL}.postfilter_summary.txt"
LOG_FILE="$LOG_DIR/${SAMPLE_FULL}.postfilter.log"

# ===================================================
# Create output directories (if they still don't exist)
# ===================================================
mkdir -p "$QC_DIR/raw"
mkdir -p "$QC_DIR/trimmed"
mkdir -p "$QC_DIR/md_flagstat"
mkdir -p "$TRIM_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$ALIGN_DIR"

# ============================================================
# Run FastQC
# ============================================================

fastqc \
  --threads "$THREADS" \
  --outdir "$QC_DIR/raw" \
  "$RAW_FASTQ_DIR"/*.fastq.gz

echo "FastQC completed successfully."

# ============================================================
# Run MultiQC
# ============================================================

multiqc "$QC_DIR/raw" --outdir "$QC_DIR/raw"

echo "MultiQC completed successfully."

# ============================================================
# Cutadapt trimming/filtering + full QC report
# ============================================================

cutadapt \
  -u 5 -u -5 \
  -U 5 -U -5 \
  -q 20,20 \
  -m 30 \
  -a A{10} \
  -A A{10} \
  -j "$THREADS" \
  --report=full \
  -o "$TRIM_DIR/${SAMPLE_FULL}_R1.trimmed.fastq.gz" \
  -p "$TRIM_DIR/${SAMPLE_FULL}_R2.trimmed.fastq.gz" \
  "$RAW_FASTQ_DIR/${SAMPLE}_1.fastq.gz" \
  "$RAW_FASTQ_DIR/${SAMPLE}_2.fastq.gz" \
  > "$LOG_DIR/cutadapt_${SAMPLE_FULL}.log"

echo "Cutadapt completed successfully."

# --- Run FastQC ---
fastqc \
  --threads "$THREADS" \
  --outdir "$QC_DIR/trimmed" \
  "$TRIM_DIR"/*.fastq.gz

echo "Post trimming FastQC completed successfully."

# --- Run MultiQC ---
multiqc "$QC_DIR/trimmed" --outdir "$QC_DIR/trimmed"

echo "Post trimming MultiQC completed successfully."


# =========================
# Step 0: Check / build BWA index
# =========================

echo "Checking BWA index..."

if [[ ! -f "${REF_FASTA}.64.bwt" && ! -f "${REF_FASTA}.bwt" ]]; then
  echo "BWA index not found for reference genome."
  echo "Building BWA index for reference genome..."
  bwa index "$REF_FASTA" 2> "$LOG_DIR/bwa_index.log"
  echo "BWA indexing completed."
else
  echo "BWA index found. Skipping indexing of reference genome."
fi

# =========================
# Step 1: Alignment: BWA-MEM with Read Groups → SAM
# =========================

echo "Running BWA-MEM alignment..."

bwa mem \
  -t "$THREADS" \
  -R "@RG\tID:${RG_ID}\tSM:${RG_SM}\tLB:${RG_LB}\tPL:${RG_PL}\tPU:${RG_PU}" \
  "$REF_FASTA" \
  "$TRIM_DIR/${SAMPLE_FULL}_R1.trimmed.fastq.gz" \
  "$TRIM_DIR/${SAMPLE_FULL}_R2.trimmed.fastq.gz" \
  > "$ALIGN_DIR/${SAMPLE_FULL}.sam" \
  2> "$LOG_DIR/bwa_mem.log"

echo "Alignment completed."

# =========================
# Step 2: Convert SAM → BAM
# =========================

if [[ ! -s "$ALIGN_DIR/${SAMPLE_FULL}.sam" ]]; then
  echo "ERROR: SAM file not created!" >&2
  exit 1
fi

echo "Converting SAM to BAM..."

samtools view -@ "$THREADS" -bS \
  "$ALIGN_DIR/${SAMPLE_FULL}.sam" \
  > "$ALIGN_DIR/${SAMPLE_FULL}.bam"

# =========================
# Step 3: Sort BAM (coordinate sort)
# =========================

if [[ ! -s "$ALIGN_DIR/${SAMPLE_FULL}.bam" ]]; then
  echo "ERROR: BAM file not created!" >&2
  exit 1
fi

echo "Sorting BAM..."

samtools sort -@ "$THREADS" \
  -o "$ALIGN_DIR/${SAMPLE_FULL}.sorted.bam" \
  "$ALIGN_DIR/${SAMPLE_FULL}.bam"

rm "$ALIGN_DIR/${SAMPLE_FULL}.sam" "$ALIGN_DIR/${SAMPLE_FULL}.bam"

echo "Sorting completed."

# =========================
# Step 4: Mark duplicates (AMPLICON-AWARE TAGGING OF DUPLICATES)
# =========================

if [[ ! -s "$ALIGN_DIR/${SAMPLE_FULL}.sorted.bam" ]]; then
  echo "ERROR: Sorted BAM not created!" >&2
  exit 1
fi

echo "Marking duplicates..."

picard MarkDuplicates \
  INPUT="$ALIGN_DIR/${SAMPLE_FULL}.sorted.bam" \
  OUTPUT="$ALIGN_DIR/${SAMPLE_FULL}.sorted.markdup.bam" \
  METRICS_FILE="$ALIGN_DIR/${SAMPLE_FULL}.markdup.metrics.txt" \
  CREATE_INDEX=false \
  REMOVE_DUPLICATES=false \
  TAG_DUPLICATE_SET_MEMBERS=true \
  VALIDATION_STRINGENCY=SILENT \
  2> "$LOG_DIR/markduplicates.log"

echo "Duplicate marking completed."

# =========================
# Step 5: Add MD and NM tags (GATK robustness) + Index final BAM (REQUIRED for GATK)
# =========================

if [[ ! -s "$ALIGN_DIR/${SAMPLE_FULL}.sorted.markdup.bam" ]]; then
  echo "ERROR: MarkDuplicates failed!" >&2
  exit 1
fi

echo "Adding MD tags..."

samtools calmd -b \
  "$ALIGN_DIR/${SAMPLE_FULL}.sorted.markdup.bam" \
  "$REF_FASTA" \
  > "$FINAL_BAM"

if [[ ! -s "$FINAL_BAM" ]]; then
  echo "ERROR: Final BAM not created by samtools calmd!" >&2
  exit 1
fi

echo "Indexing final BAM..."

samtools index "$FINAL_BAM"

echo "BAM indexing completed."

# =========================
# Step 6: Alignment statistics
# =========================

echo "Generating alignment statistics..."

samtools flagstat \
  "$FINAL_BAM" \
  > "$LOG_DIR/${SAMPLE_FULL}.flagstat.txt"

# =========================
# Step 7: MultiQC (duplicates + alignment metrics)
# =========================

echo "Running MultiQC..."

multiqc \
  "$ALIGN_DIR" \
  "$LOG_DIR" \
  --outdir "$QC_DIR/md_flagstat"

echo "MultiQC of MarkDuplicates & flagstat & Cutadapt done!"

# =========================
# Step 8: Cleanup intermediate files (Optional)
# =========================

echo "Cleaning up intermediate files..."

rm -f \
  "$ALIGN_DIR/${SAMPLE_FULL}.sorted.bam" \
  "$ALIGN_DIR/${SAMPLE_FULL}.sorted.markdup.bam"

echo "Cleanup completed."


echo "Alignment and BAM preprocessing completed successfully."

# ========================
# Time elapsed
# ========================

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "----------------------------------------"
echo "Alignment + BAM preprocessing completed."
echo "Total runtime: ${ELAPSED} seconds"
echo "Total runtime: $((ELAPSED / 60)) minutes"
echo "----------------------------------------"
echo " "
echo "================ PREPROCESSING DONE ================"
echo " "

# ========================
# Somatic variant calling
# ========================

START_TIME_MUTECT2=$(date +%s)

echo "Starting somatic variant calling with Mutect2..."

# ===================================================
# Create output directories (if they still don't exist)
# ===================================================

mkdir -p "$VARIANT_DIR"

# ============================================================
# Sanity checks (fail early, fail clearly)
# ============================================================

echo "Running sanity checks..."

for file in \
  "$REF_FASTA" \
  "$REF_FASTA_DICT" \
  "$REF_FASTA.fai" \
  "$INPUT_MD_BAM" \
  "$INPUT_MD_BAM.bai" \
  "$INTERVALS" \
  "$PON" \
  "$PON.tbi" \
  "$GNOMAD" \
  "$GNOMAD.tbi"
do
  [[ -s "$file" ]] || { echo "ERROR: Missing required file: $file"; exit 1; }
done

echo "All required input files found."

# ============================================================
# Step 1: Mutect2 (tumor-only panel of 7 genes)
# ============================================================

echo "Running Mutect2..."

gatk --java-options "$JAVA_MEM" Mutect2 \
  -R "$REF_FASTA" \
  -I "$INPUT_MD_BAM" \
  --tumor-sample "$TUMOR_SM" \
  --panel-of-normals "$PON" \
  --germline-resource "$GNOMAD" \
  -L "$INTERVALS" \
  --af-of-alleles-not-in-resource 0.0000025 \
  --f1r2-tar-gz "$VARIANT_DIR/${SAMPLE_FULL}.f1r2.tar.gz" \
  -O "$VARIANT_DIR/${SAMPLE_FULL}.unfiltered.vcf.gz" \
  > "$LOG_DIR/mutect2.stdout.log" \
  2> "$LOG_DIR/mutect2.stderr.log"


echo "Mutect2 completed successfully."
echo "Unfiltered VCF written to:"
echo "  $VARIANT_DIR/${SAMPLE_FULL}.unfiltered.vcf.gz"

# ========================
# Time elapsed
# ========================

END_TIME_2=$(date +%s)
ELAPSED_2=$((END_TIME_2 - START_TIME_MUTECT2))

echo "----------------------------------------"
echo "Somatic variant calling completed."
echo "Total runtime: ${ELAPSED_2} seconds"
echo "Total runtime: $((ELAPSED_2 / 60)) minutes"
echo "----------------------------------------"
echo " "

echo "Learning read orientation model starts after sanity checks..."

# ============================================================
# Sanity checks
# ============================================================

echo "Running sanity checks..."

# Sanity checks
for file in \
  "$F1R2_TAR" \
  "$INPUT_MD_BAM" \
  "$INPUT_MD_BAM.bai" \
  "$REF_FASTA" \
  "$REF_FASTA.fai" \
  "$REF_FASTA_DICT" \
  "$GNOMAD" \
  "$GNOMAD.tbi"
do
  [[ -s "$file" ]] || { echo "ERROR: Missing required file: $file"; exit 1; }
done

echo "All required files found."

# ============================================================
# LearnReadOrientationModel
# ============================================================

echo "Learning read orientation model..."

gatk --java-options "$JAVA_MEM" LearnReadOrientationModel \
  -I "$F1R2_TAR" \
  -O "$ORIENTATION_MODEL" \
  2> "$LOG_DIR/learn_read_orientation_model.log"

echo "LearnReadOrientationModel completed successfully."
echo "Orientation model written to:"
echo "  $ORIENTATION_MODEL"

# If $ORIENTATION_MODEL output is empty -> "Error"
[[ -s "$ORIENTATION_MODEL" ]] || { echo "Orientation model failed"; exit 1; }

# ============================================================
# GetPileupSummaries
# ============================================================
echo " "
echo "Starting GetPileupSummaries..."

gatk --java-options "$JAVA_MEM" GetPileupSummaries \
  -R "$REF_FASTA" \
  -I "$INPUT_MD_BAM" \
  -V "$GNOMAD" \
  -L "$INTERVALS" \
  -O "$PILEUP_TABLE" \
  2> "$LOG_DIR/get_pileup_summaries.log"

echo "GetPileupSummaries completed."

# To get to know the amount of informative SNPs. If number of sites < 10–20, contamination estimate is statistically weak.
echo " "

SITE_COUNT=$(grep -v "^#" "$PILEUP_TABLE" | wc -l)

if [[ "$SITE_COUNT" -lt 10 ]]; then
  echo "WARNING: Very few informative SNPs ($SITE_COUNT). Contamination estimate may be unreliable."
fi


echo "CalculateContamination will start after sanity checks"
# ============================================================
# Sanity checks
# ============================================================

echo "Running sanity checks..."

for file in "$PILEUP_TABLE"
do
  [[ -s "$file" ]] || { echo "ERROR: Missing $file"; exit 1; }
done

echo "All required files found."

# ============================================================
# CalculateContamination
# ============================================================
echo "Starting CalculateContamination..."

gatk --java-options "$JAVA_MEM" CalculateContamination \
  -I "$PILEUP_TABLE" \
  -O "$CONTAM_TABLE" \
  2> "$LOG_DIR/calculate_contamination.log"

echo "CalculateContamination completed."
echo "Contamination table written to:"
echo "  $CONTAM_TABLE"

# If $CONTAM_TABLE output is empty -> "Error"
[[ -s "$CONTAM_TABLE" ]] || { echo "Contamination estimation failed"; exit 1; }

echo " "
echo "Variant filtering will start after sanity checks"
# ============================================================
# Sanity checks
# ============================================================

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


echo "All required files found."

# ============================================================
# FilterMutectCalls
# ============================================================

echo "Filtering Mutect2 calls..."

gatk --java-options "$JAVA_MEM" FilterMutectCalls \
  -R "$REF_FASTA" \
  -V "$UNFILTERED_VCF" \
  --contamination-table "$CONTAM_TABLE" \
  --orientation-bias-artifact-priors "$ORIENTATION_MODEL" \
  -O "$FILTERED_VCF" \
  2> "$LOG_DIR/filter_mutect_calls.log"

echo "FilterMutectCalls completed successfully."
echo "Final filtered VCF:"
echo "  $FILTERED_VCF"

# Check if $FILTERED_VCF output exists and is not empty -> otherwise "Error"
[[ -s "$FILTERED_VCF" ]] || { echo "FilterMutectCalls failed"; exit 1; }

# Find variants with 'PASS' in $FILTERED_VCF file in Terminal like this:
# "bcftools view -H -f PASS SRR30536566_full.filtered.vcf.gz"

# ============================================================
# Extract PASS variants only
# ============================================================

echo "Extracting PASS variants..."

bcftools view -f PASS "$FILTERED_VCF" -Oz -o "$PASS_VCF"       # "-f PASS" only shows variants where the FILTER column == PASS. It will return only high-confidence calls.
bcftools index -t "$PASS_VCF"                                  # "-t" when indexing as ".tbi", otherwise ".csi" by default.

# Count PASS variants
PASS_COUNT=$(bcftools view -H "$PASS_VCF" | wc -l)             # "-H" hides header. It reads a file that already contains only PASS variants (from the '-f PASS code')
echo "Number of PASS variants: $PASS_COUNT"

echo "PASS-only VCF written to:"
echo "  $PASS_VCF"



echo " "
echo "Post variant filtering will start after sanity checks"

# ============================================================
# Thresholds (amplicon tumor-only)
# ============================================================

MIN_DP=200          # total depth
MIN_AD_ALT=10       # ALT read count
MIN_VAF=0.02        # 2%

# ============================================================
# Logging
# ============================================================

exec > >(tee -a "$LOG_FILE") 2>&1

echo "Starting post-filtering"
echo "Sample: $SAMPLE_FULL"

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

echo "Applying post-filter thresholds:"
echo "  DP >= $MIN_DP"
echo "  ALT reads (AD[1]) >= $MIN_AD_ALT"
echo "  VAF >= $MIN_VAF"

bcftools view -f PASS "$FILTERED_VCF" | \
bcftools filter \
  -i "FORMAT/DP >= ${MIN_DP} && FORMAT/AD[0:1] >= ${MIN_AD_ALT} && FORMAT/AF >= ${MIN_VAF}" \
  -Oz -o "$POSTFILTER_VCF"

# ============================================================
# Index the post-filtered VCF (required for IGV) → .tbi
# ============================================================

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
  echo "Sample: $SAMPLE_FULL"
  echo ""
  echo "Library type: Amplicon (PCR)"
  echo "Sequencing: Tumor-only"
  echo ""
  echo "Thresholds:"
  echo "  DP >= $MIN_DP"
  echo "  ALT reads >= $MIN_AD_ALT"
  echo "  VAF >= $MIN_VAF"
  echo ""
  echo "Variants retained: $N_VARIANTS"
} > "$SUMMARY_TXT"

echo " "
echo "[$(date)] $SAMPLE_FULL post-filtering completed successfully"
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
