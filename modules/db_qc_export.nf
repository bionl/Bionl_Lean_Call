// modules/db_qc_export.nf
// ───────────────────────────────────────────────────────────────────────────
// DB QC Export — per-sample QC metrics for downstream ingestion decisions.
//
// Runs samtools flagstat, mosdepth, and bcftools stats in parallel per
// sample, then aggregates results into a structured JSON + TSV.
//
// This module is advisory only — it never blocks or filters samples.
// ───────────────────────────────────────────────────────────────────────────
nextflow.enable.dsl=2

params.outdir    = params.outdir    ?: "results"
params.scriptdir = params.scriptdir ?: "${workflow.projectDir}/scripts"

/********************  PROCESSES  ********************/

process DB_QC_FLAGSTAT {
    tag { "${meta.sample}" }

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("${meta.sample}.flagstat.txt")

    script:
    def sample = meta.sample
    """
    samtools flagstat ${bam} > ${sample}.flagstat.txt
    """
}

process DB_QC_MOSDEPTH {
    tag { "${meta.sample}" }

    input:
        tuple val(meta), path(bam), path(bai)
        path bed

    output:
        tuple val(meta),
            path("${meta.sample}.mosdepth.summary.txt"),
            path("${meta.sample}.thresholds.bed.gz"),
            path("${meta.sample}.mosdepth.global.dist.txt")

    script:
    def sample = meta.sample
    def has_bed = bed.name != 'NO_FILE'
    def bed_args = has_bed ? "--by ${bed} --thresholds 10,20,30,50" : ""
    """
    ln -sf ${bam} ${sample}.bam
    ln -sf ${bai} ${sample}.bam.bai

    mosdepth \\
        --no-per-base \\
        ${bed_args} \\
        --fast-mode \\
        ${sample} ${sample}.bam

    # Ensure thresholds file exists even for WGS (no --by)
    if [ ! -f ${sample}.thresholds.bed.gz ]; then
        echo -n | gzip > ${sample}.thresholds.bed.gz
    fi
    """
}

process DB_QC_VCFSTATS {
    tag { "${meta.sample}" }

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta),
            path("${meta.sample}.bcftools_stats.txt"),
            path("${meta.sample}.callers.txt")

    script:
    def sample = meta.sample
    """
    if [ ! -f ${vcf}.tbi ] && [ ! -f ${vcf}.csi ]; then
        bcftools index -t ${vcf} 2>/dev/null || true
    fi

    bcftools stats ${vcf} > ${sample}.bcftools_stats.txt

    if bcftools view -h ${vcf} | grep -q 'ID=CALLERS'; then
        bcftools query -f '%INFO/CALLERS\\n' ${vcf} > ${sample}.callers.txt
    else
        touch ${sample}.callers.txt
    fi
    """
}

process DB_QC_EXPORT_JSON {
    tag { "${meta.sample}" }
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
        tuple val(meta),
            path(flagstat),
            path(mosdepth_summary),
            path(mosdepth_thresholds),
            path(mosdepth_global_dist),
            path(bcftools_stats),
            path(callers)
        each path(script)

    output:
        tuple val(meta), path("${meta.sample}.qc.json"), path("${meta.sample}.qc.tsv")

    script:
    def sample = meta.sample
    def assay  = meta.assay ?: "NA"
    """
    python ${script} \\
        --flagstat ${flagstat} \\
        --mosdepth-summary ${mosdepth_summary} \\
        --mosdepth-thresholds ${mosdepth_thresholds} \\
        --mosdepth-global-dist ${mosdepth_global_dist} \\
        --bcftools-stats ${bcftools_stats} \\
        --callers ${callers} \\
        --sample ${sample} \\
        --assay ${assay} \\
        --outdir .
    """
}

process DB_QC_AGGREGATE {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
        path(tsvs)

    output:
        path("qc_summary.tsv")

    script:
    """
    head -1 \$(ls *.qc.tsv | head -1) > qc_summary.tsv
    tail -q -n +2 *.qc.tsv >> qc_summary.tsv
    """
}

/********************  WORKFLOW  ********************/

workflow DB_QC_EXPORT {
    take:
        vcf_ch   // tuple(meta, vcf)
        bam_ch   // tuple(meta, bam, bai)
        bed_ch   // value channel with BED file

    main:
        qc_script_ch = Channel
            .fromPath("${params.scriptdir}/compute_qc_metrics.py")
            .first()

        // Three independent per-sample analyses
        DB_QC_FLAGSTAT(bam_ch)
        DB_QC_MOSDEPTH(bam_ch, bed_ch)
        DB_QC_VCFSTATS(vcf_ch)

        // Join all outputs by meta key → single tuple per sample
        qc_input_ch = DB_QC_FLAGSTAT.out                          // (meta, flagstat)
            .join(DB_QC_MOSDEPTH.out)                             // + (summary, thresholds, global_dist)
            .join(DB_QC_VCFSTATS.out)                             // + (bcfstats, callers)

        DB_QC_EXPORT_JSON(qc_input_ch, qc_script_ch)

        // Aggregate all per-sample TSVs into one summary
        all_tsvs = DB_QC_EXPORT_JSON.out
            .map { meta, json, tsv -> tsv }
            .collect()

        DB_QC_AGGREGATE(all_tsvs)

    emit:
        qc_json    = DB_QC_EXPORT_JSON.out   // tuple(meta, json, tsv) per sample
        qc_summary = DB_QC_AGGREGATE.out     // single aggregated TSV
}
