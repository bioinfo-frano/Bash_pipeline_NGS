nextflow.enable.dsl=2

// ------------------------------------------------------------
// PARAMETERS
// ------------------------------------------------------------

// Replace "/path/to/your/Genomics_cancer" with the actual full path to your "Genomics_cancer" folder.
// Example: 
//   on macOS: /Users/john/Desktop/Projects/Genomics_cancer
//   on Linux: /home/john/Projects/Genomics_cancer

base_dir = "/path/to/your/Genomics_cancer"

params.fastq = "$base_dir/data/SRR30536566/raw_fastq/*_{1,2}.fastq.gz"

params.qc_report = "$base_dir/data/SRR30536566_full_nf/qc"
params.trimming  = "$base_dir/data/SRR30536566_full_nf/trimmed"
params.logs      = "$base_dir/data/SRR30536566_full_nf/logs"
params.aligned   = "$base_dir/data/SRR30536566_full_nf/aligned"
params.variants  = "$base_dir/data/SRR30536566_full_nf/variants"

params.reference_fasta = "$base_dir/reference/GRCh38/fasta/Homo_sapiens_assembly38.fasta"
params.intervals = "$base_dir/reference/GRCh38/intervals/crc_panel_7genes_sorted.hg38.bed"
params.pon       = "$base_dir/reference/GRCh38/somatic_resources/1000g_pon.hg38.vcf.gz"
params.gnomad    = "$base_dir/reference/GRCh38/somatic_resources/af-only-gnomad.hg38.vcf.gz"

// ------------------------------------------------------------
// CHANNEL: read pairs
// ------------------------------------------------------------
fastq_ch = Channel.fromFilePairs(params.fastq, flat:true)
          .map { sample_id, r1, r2 -> tuple(sample_id, [r1, r2]) }

// ------------------------------------------------------------
// PREPROCESSING PIPELINE
// ------------------------------------------------------------
process FASTQC {
    conda "DNA2"
    tag "$sample_id"
    cpus 4
    publishDir "${params.qc_report}/raw", mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html")

    script:
    """
    fastqc -t ${task.cpus} ${reads.join(' ')}
    """
}

process MULTIQC_RAW {
    conda "DNA2"
    cpus 2
    publishDir "${params.qc_report}/raw", mode:'copy'

    input:
    path fastqc_zips

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}

process CUTADAPT {
    conda "DNA2"
    tag "$sample_id"
    cpus 4

    publishDir "${params.trimming}", mode:'copy', pattern:"*trimmed.fastq.gz"
    publishDir "${params.logs}", mode:'copy', pattern:"cutadapt*.log"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("*R1.trimmed.fastq.gz"),
          path("*R2.trimmed.fastq.gz"),
          path("cutadapt_*.log")

    script:
    """
    THREADS=${task.cpus}
    R1_TRIM="${sample_id}_full_nf_R1.trimmed.fastq.gz"
    R2_TRIM="${sample_id}_full_nf_R2.trimmed.fastq.gz"
    LOG_FILE="cutadapt_${sample_id}_full_nf.log"

    cutadapt -u 5 -u -5 -U 5 -U -5 -q 20,20 -m 30 -a A{10} -A A{10} -j \$THREADS \
      --report=full -o \$R1_TRIM -p \$R2_TRIM ${reads[0]} ${reads[1]} \
      > \$LOG_FILE 2>&1
    """
}

process FASTQC_TRIMMED {
    conda "DNA2"
    tag "$sample_id"
    cpus 4
    publishDir "${params.qc_report}/trimmed", mode:'copy'

    input:
    tuple val(sample_id), path(r1_trim), path(r2_trim)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html")

    script:
    """
    fastqc -t ${task.cpus} ${r1_trim} ${r2_trim}
    """
}

process MULTIQC_TRIMMED {
    conda "DNA2"
    cpus 2
    publishDir "${params.qc_report}/trimmed", mode:'copy'

    input:
    path fastqc_zips

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}

// ------------------------------------------------------------
// ALIGNMENT
// ------------------------------------------------------------
process ALIGNMENT {
    conda "DNA2"
    tag "$sample_id"
    cpus 4
    memory '8 GB'

    publishDir "${params.aligned}", mode:'copy', pattern:"*.bam"
    publishDir "${params.aligned}", mode:'copy', pattern:"*.bai"
    publishDir "${params.aligned}", mode:'copy', pattern:"*.metrics.txt"
    publishDir "${params.logs}", mode:'copy', pattern:"*.log"
    publishDir "${params.logs}", mode:'copy', pattern:"*.flagstat.txt"

    input:
    tuple val(sample_id), path(r1_trim), path(r2_trim), path(cutadapt_log)

    output:
    tuple val(sample_id),
          path("${sample_id}_full_nf.sorted.markdup.md.bam"),
          path("${sample_id}_full_nf.sorted.markdup.md.bam.bai"),
          path("${sample_id}_full_nf.markdup.metrics.txt"),
          path("${sample_id}_full_nf.flagstat.txt"),
          path(cutadapt_log),
          path("bwa_mem.log"),
          path("markduplicates.log")

    script:
    """
    set -euo pipefail
    THREADS=${task.cpus}
    SAMPLE="${sample_id}"
    SAMPLE_FULL="\${SAMPLE}_full_nf"
    REF_FASTA="${params.reference_fasta}"

    bwa mem -t 2 -R "@RG\\tID:\${SAMPLE}\\tSM:DMBEL-EIDR-071\\tLB:AMPLICON\\tPL:ILLUMINA\\tPU:HiSeq4000" \
        "\${REF_FASTA}" ${r1_trim} ${r2_trim} > "\${SAMPLE_FULL}.sam" 2> bwa_mem.log

    samtools view -@ \$THREADS -b "\${SAMPLE_FULL}.sam" > "\${SAMPLE_FULL}.bam"
    samtools sort -@ \$THREADS -m 2G -o "\${SAMPLE_FULL}.sorted.bam" "\${SAMPLE_FULL}.bam"
    rm "\${SAMPLE_FULL}.sam" "\${SAMPLE_FULL}.bam"

    picard MarkDuplicates \
        INPUT="\${SAMPLE_FULL}.sorted.bam" \
        OUTPUT="\${SAMPLE_FULL}.sorted.markdup.bam" \
        METRICS_FILE="\${SAMPLE_FULL}.markdup.metrics.txt" \
        CREATE_INDEX=false REMOVE_DUPLICATES=false TAG_DUPLICATE_SET_MEMBERS=true \
        VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true 2>&1 | tee markduplicates.log

    samtools calmd -b "\${SAMPLE_FULL}.sorted.markdup.bam" "\${REF_FASTA}" > "\${SAMPLE_FULL}.sorted.markdup.md.bam"
    samtools index "\${SAMPLE_FULL}.sorted.markdup.md.bam"
    samtools flagstat "\${SAMPLE_FULL}.sorted.markdup.md.bam" > "\${SAMPLE_FULL}.flagstat.txt"
    rm "\${SAMPLE_FULL}.sorted.bam" "\${SAMPLE_FULL}.sorted.markdup.bam"
    """
}

process MULTIQC_ALIGNMENT {
    conda "DNA2"
    cpus 2
    publishDir "${params.qc_report}/md_flagstat", mode:'copy'

    input:
    tuple path(metrics), path(flagstat), path(cutadapt_log), path(bwa_log), path(markdup_log)

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}

// ------------------------------------------------------------
// VARIANT CALLING
// ------------------------------------------------------------
process MUTECT2 {
    conda "DNA2"
    cpus 4
    memory '8 GB'
    publishDir "${params.variants}", mode:'copy', pattern:"*.vcf.gz*"
    publishDir "${params.variants}", mode:'copy', pattern:"*.f1r2.tar.gz"
    publishDir "${params.logs}", mode:'copy', pattern:"mutect2*.log"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
          path("${sample_id}_full_nf.unfiltered.vcf.gz"),
          path("${sample_id}_full_nf.unfiltered.vcf.gz.tbi"),
          path("${sample_id}_full_nf.unfiltered.vcf.gz.stats"),
          path("${sample_id}_full_nf.f1r2.tar.gz"),
          path("mutect2.stdout.log"),
          path("mutect2.stderr.log")

    script:
    """
    SAMPLE_FULL="${sample_id}_full_nf"
    gatk --java-options "-Xmx6g" Mutect2 \
      -R ${params.reference_fasta} \
      -I ${bam} \
      --tumor-sample DMBEL-EIDR-071 \
      --panel-of-normals ${params.pon} \
      --germline-resource ${params.gnomad} \
      -L ${params.intervals} \
      --af-of-alleles-not-in-resource 0.0000025 \
      --f1r2-tar-gz \${SAMPLE_FULL}.f1r2.tar.gz \
      -O \${SAMPLE_FULL}.unfiltered.vcf.gz \
      > mutect2.stdout.log 2> mutect2.stderr.log
    """
}

process LEARN_READ_ORIENTATION {
    conda "DNA2"
    publishDir "${params.variants}", mode:'copy', pattern:"*.read-orientation-model.tar.gz"
    publishDir "${params.logs}", mode:'copy', pattern:"*.log"      // ← changed to wildcard

    input:
    tuple val(sample_id), path(f1r2_tar)

    output:
    tuple val(sample_id),
          path("${sample_id}_full_nf.read-orientation-model.tar.gz"),
          path("learn_read_orientation_model.log")

    script:
    """
    SAMPLE_FULL="${sample_id}_full_nf"
    gatk LearnReadOrientationModel \\
      -I ${f1r2_tar} \\
      -O \${SAMPLE_FULL}.read-orientation-model.tar.gz \\
      > learn_read_orientation_model.log 2>&1
    """
}

process GET_PILEUP_SUMMARIES {
    conda "DNA2"
    publishDir "${params.variants}", mode:'copy', pattern:"*.pileups.table"
    publishDir "${params.logs}", mode:'copy', pattern:"*.log"      // ← changed to wildcard

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
          path("${sample_id}_full_nf.pileups.table"),
          path("get_pileup_summaries.log")

    script:
    """
    SAMPLE_FULL="${sample_id}_full_nf"
    gatk GetPileupSummaries \\
      -R ${params.reference_fasta} \\
      -I ${bam} \\
      -V ${params.gnomad} \\
      -L ${params.intervals} \\
      -O \${SAMPLE_FULL}.pileups.table \\
      > get_pileup_summaries.log 2>&1
    """
}

process CALCULATE_CONTAMINATION {
    conda "DNA2"
    publishDir "${params.variants}", mode:'copy', pattern:"*.contamination.table"
    publishDir "${params.logs}", mode:'copy', pattern:"*.log"      // ← changed to wildcard

    input:
    tuple val(sample_id), path(pileups_table)

    output:
    tuple val(sample_id),
          path("${sample_id}_full_nf.contamination.table"),
          path("calculate_contamination.log")

    script:
    """
    SAMPLE_FULL="${sample_id}_full_nf"
    gatk CalculateContamination \\
      -I ${pileups_table} \\
      -O \${SAMPLE_FULL}.contamination.table \\
      > calculate_contamination.log 2>&1
    """
}

process FILTER_MUTECT_CALLS {
    conda "DNA2"
    publishDir "${params.variants}", mode:'copy', pattern:"*.filtered.vcf.gz*"
    publishDir "${params.logs}", mode:'copy', pattern:"*.log"      // ← changed to wildcard

    input:
    tuple val(sample_id), path(unfiltered_vcf), path(unfiltered_tbi), path(unfiltered_stats),
          path(orientation_model), path(contamination_table)

    output:
    tuple val(sample_id),
          path("${sample_id}_full_nf.filtered.vcf.gz"),
          path("${sample_id}_full_nf.filtered.vcf.gz.tbi"),
          path("${sample_id}_full_nf.filtered.vcf.gz.filteringStats.tsv"),
          path("filter_mutect_calls.log")

    script:
    """
    SAMPLE_FULL="${sample_id}_full_nf"
    gatk FilterMutectCalls \\
      -R ${params.reference_fasta} \\
      -V ${unfiltered_vcf} \\
      --contamination-table ${contamination_table} \\
      --orientation-bias-artifact-priors ${orientation_model} \\
      -O \${SAMPLE_FULL}.filtered.vcf.gz \\
      > filter_mutect_calls.log 2>&1
    """
}

// --- POSTFILTER_VARIANTS ---
process POSTFILTER_VARIANTS {
    conda "DNA2"
    publishDir "${params.variants}", mode:'copy', pattern:"*.postfiltered.vcf.gz*"
    publishDir "${params.variants}", mode:'copy', pattern:"*.filtered.PASS.vcf.gz*"
    publishDir "${params.variants}", mode:'copy', pattern:"*.postfilter_summary.txt"
    publishDir "${params.logs}", mode:'copy', pattern:"*.log"      // ← already wildcard, but kept

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_tbi)

    output:
    path("${sample_id}_full_nf.postfiltered.vcf.gz")
    path("${sample_id}_full_nf.postfiltered.vcf.gz.tbi")
    path("${sample_id}_full_nf.filtered.PASS.vcf.gz")
    path("${sample_id}_full_nf.filtered.PASS.vcf.gz.tbi")
    path("${sample_id}_full_nf.postfilter_summary.txt")
    path("${sample_id}_full_nf.postfilter.log")

    script:
    """
    SAMPLE_FULL="${sample_id}_full_nf"
    MIN_DP=200
    MIN_AD_ALT=10
    MIN_VAF=0.02

    # Start log
    {
        echo "Starting post-filtering"
        echo "Sample: \$SAMPLE_FULL"
        echo "Applying post-filter thresholds:"
        echo "  DP >= \$MIN_DP"
        echo "  ALT reads (AD[1]) >= \$MIN_AD_ALT"
        echo "  VAF >= \$MIN_VAF"
    } > \${SAMPLE_FULL}.postfilter.log

    # Extract PASS variants
    bcftools view -f PASS ${filtered_vcf} -Oz -o \${SAMPLE_FULL}.filtered.PASS.vcf.gz
    bcftools index -t \${SAMPLE_FULL}.filtered.PASS.vcf.gz

    # Apply hard filters
    bcftools view -f PASS ${filtered_vcf} | \\
      bcftools filter -i "FORMAT/DP >= \${MIN_DP} && FORMAT/AD[0:1] >= \${MIN_AD_ALT} && FORMAT/AF >= \${MIN_VAF}" \\
      -Oz -o \${SAMPLE_FULL}.postfiltered.vcf.gz
    bcftools index -t \${SAMPLE_FULL}.postfiltered.vcf.gz

    # Count variants
    N_VARIANTS=\$(bcftools view -H \${SAMPLE_FULL}.postfiltered.vcf.gz | wc -l)

    # Continue log
    {
        echo "Indexing post-filtered VCF"
        echo "Variants retained after post-filtering: \$N_VARIANTS"
        echo ""
        echo "\$(date) \$SAMPLE_FULL post-filtering completed successfully"
    } >> \${SAMPLE_FULL}.postfilter.log

    # Also generate bcftools stats summary file (optional)
    bcftools stats \${SAMPLE_FULL}.postfiltered.vcf.gz > \${SAMPLE_FULL}.postfilter_summary.txt
    """
}

// ------------------------------------------------------------
// WORKFLOW
// ------------------------------------------------------------
workflow {

    // ---- Raw FASTQC + MultiQC ----
    fastqc_raw = FASTQC(fastq_ch)
    fastqc_raw_zips = fastqc_raw.map { sample_id, zips, htmls -> zips }.flatten()
    MULTIQC_RAW( fastqc_raw_zips.collect() )

    // ---- Trimming ----
    trimmed = CUTADAPT(fastq_ch)

    // ---- Trimmed FASTQC + MultiQC ----
    fastqc_trim = FASTQC_TRIMMED( trimmed.map { [it[0], it[1], it[2]] } )
    fastqc_trim_zips = fastqc_trim.map { sample_id, zips, htmls -> zips }.flatten()
    MULTIQC_TRIMMED( fastqc_trim_zips.collect() )

    // ---- Alignment ----
    aligned = ALIGNMENT(trimmed)

    // ---- MultiQC for alignment ----
    alignment_for_multiqc = aligned.map { sample_id, bam, bai, metrics, flagstat, cutadapt_log, bwa_log, markdup_log ->
        tuple(metrics, flagstat, cutadapt_log, bwa_log, markdup_log)
    }
    MULTIQC_ALIGNMENT(alignment_for_multiqc)

    // --- Variant calling ---
    // Channel with BAM only
    bam_ch = aligned.map { sample_id, bam, bai, metrics, flagstat, cutadapt_log, bwa_log, markdup_log ->
        tuple(sample_id, bam, bai)
    }

    // Mutect2
    mutect_out = MUTECT2(bam_ch)

    // Split mutect_out into separate channels for clarity
    mutect_for_orient = mutect_out.map { sample_id, vcf, tbi, stats, f1r2, stdout, stderr ->
        tuple(sample_id, f1r2)
    }
    mutect_for_filter = mutect_out.map { sample_id, vcf, tbi, stats, f1r2, stdout, stderr ->
        tuple(sample_id, vcf, tbi, stats)
    }

    // Orientation model (strip log before joining)
    orient_out_raw = LEARN_READ_ORIENTATION(mutect_for_orient)
    orient_for_join = orient_out_raw.map { sample_id, model, log -> tuple(sample_id, model) }

    // Pileups (strip log before passing to contamination)
    pileups_raw = GET_PILEUP_SUMMARIES(bam_ch)
    pileups_for_contam = pileups_raw.map { sample_id, table, log -> tuple(sample_id, table) }

    // Contamination (strip log before joining)
    contam_raw = CALCULATE_CONTAMINATION(pileups_for_contam)
    contam_for_join = contam_raw.map { sample_id, table, log -> tuple(sample_id, table) }

    // Combine inputs for FilterMutectCalls
    filter_input = mutect_for_filter
        .join(orient_for_join, by: 0)
        .join(contam_for_join, by: 0)
        .map { sample_id, vcf, tbi, stats, orient_model, contam_table ->
            tuple(sample_id, vcf, tbi, stats, orient_model, contam_table)
        }

    filter_out = FILTER_MUTECT_CALLS(filter_input)

    // Final post‑filtering
    POSTFILTER_VARIANTS( filter_out.map { sample_id, filtered_vcf, filtered_tbi, stats, filter_log ->
        tuple(sample_id, filtered_vcf, filtered_tbi)
    } )
}