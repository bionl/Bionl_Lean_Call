nextflow.enable.dsl=2

// ═══════════════════════════════════════════════════════════════════════════
// IMPORTS
// ═══════════════════════════════════════════════════════════════════════════

include { PIPELINE_INITIALISATION; PIPELINE_COMPLETION } \
  from './external/sarek/subworkflows/local/utils_nfcore_sarek_pipeline'

include { NFCORE_SAREK } \
  from './external/sarek/main.nf'

include { POST_SAREK } \
  from './modules/vep.nf'

include { CONSENSUS_CALLING } \
  from './modules/consensus.nf'

// ═══════════════════════════════════════════════════════════════════════════
// PARAMETERS & VALIDATION
// ═══════════════════════════════════════════════════════════════════════════

// Default parameters
params.input                  = params.input ?: params.samplesheet
params.outdir                 = params.outdir ?: params.output
params.bed                    = params.bed ?: "${workflow.projectDir}/data/annotated_merged_MANE_deduped.bed"
//params.run_variant_calling    = params.run_variant_calling instanceof Boolean ? params.run_variant_calling : true
params.create_consensus       = params.create_consensus instanceof Boolean ? params.create_consensus : true
params.ref_fasta              = params.ref_fasta ?: params.vep_fasta
//params.vep_fasta              = params.vep_fasta ?: params.vep_fasta

// Validate required parameters
//if (params.run_variant_calling) {
//    if (!params.input)  error "❌ Missing --input (samplesheet CSV) when run_variant_calling=true"
//    if (!params.outdir) error "❌ Missing --outdir when run_variant_calling=true"
//    if (params.create_consensus && !params.ref_fasta) {
//        error "❌ Missing --vep_fasta when create_consensus=true"
//    }
//} else {
//    if (!params.post_samplesheet && !params.variant_calling_outdir)
//        error "❌ When run_variant_calling=false provide either --post_samplesheet or --variant_calling_outdir"
//    
//    if (params.post_samplesheet && params.variant_calling_outdir)
//        error "❌ Cannot provide both --post_samplesheet and --variant_calling_outdir. Choose one."
//}

// ═══════════════════════════════════════════════════════════════════════════
// HELPER FUNCTIONS
// ═══════════════════════════════════════════════════════════════════════════

def isGcsPath(path) {
    return path.toString().startsWith('gs://')
}

def validateBedFile() {
    def bedFile = params.bed ? file(params.bed) : null
    if (!bedFile?.exists()) {
        error "❌ BED file not found: ${params.bed}"
    }
    return bedFile
}

// ═══════════════════════════════════════════════════════════════════════════
// WORKFLOWS
// ═══════════════════════════════════════════════════════════════════════════

workflow COLLECT_VARIANT_CALLING_OUTPUTS {
    take:
        trigger    // Completion signal channel
        outdir     // Output directory to search

    main:
        def isGCS = isGcsPath(outdir)
        
        // Collect DeepVariant VCFs
        dv_vcf_ch = trigger
            .flatMap { 
                file("${outdir}/variant_calling/deepvariant/*/*.vcf.gz", checkIfExists: !isGCS)
            }
            .filter { vcf -> 
                vcf.name.endsWith('.vcf.gz') && 
                !vcf.name.contains('.g.vcf.gz') && 
                !vcf.name.endsWith('.tbi') 
            }
            .map { vcf -> tuple(vcf.parent.name, vcf) }
        
        // Collect HaplotypeCaller VCFs
        hc_vcf_ch = trigger
            .flatMap { 
                file("${outdir}/variant_calling/haplotypecaller/*/*.vcf.gz", checkIfExists: !isGCS)
            }
            .filter { vcf -> 
                vcf.name.endsWith('.vcf.gz') && 
                !vcf.name.contains('.g.vcf.gz') && 
                !vcf.name.endsWith('.tbi') 
            }
            .map { vcf -> tuple(vcf.parent.name, vcf) }

        // Collect BAMs with BAI
        bam_ch = trigger
            .flatMap { 
                file("${outdir}/preprocessing/mapped/*/*.sorted.bam", checkIfExists: !isGCS)
            }
            .map { bam -> 
                def sample = bam.parent.name
                def bamPath = bam.toString()
                def baiPath = "${bamPath}.bai"
                
                def bai
                if (isGCS) {
                    bai = file(baiPath, checkIfExists: false)
                } else {
                    bai = file(baiPath)
                    if (!bai.exists()) {
                        bai = file("${bam.parent}/${bam.baseName}.bai")
                        if (!bai.exists()) {
                            error "❌ BAM index not found for ${bam}"
                        }
                    }
                }
                
                tuple(sample, bam, bai)
            }

    emit:
        dv_vcf = dv_vcf_ch
        hc_vcf = hc_vcf_ch
        bam = bam_ch
}

workflow RUN_FROM_VARIANT_CALLING_OUTDIR {
    take:
        variant_calling_outdir
        bed_ch

    main:
        def isGCS = isGcsPath(variant_calling_outdir)
        
        log.info """
        ╔════════════════════════════════════════════════════════════╗
        ║  Using existing variant calling results                    ║
        ║  Location: ${variant_calling_outdir}
        ╚════════════════════════════════════════════════════════════╝
        """.stripIndent()

        // Collect VCFs (use consensus or single caller depending on params)
        if (params.create_consensus) {
            // Collect DeepVariant VCFs
            dv_vcf_ch = Channel
                .fromPath("${variant_calling_outdir}/variant_calling/deepvariant/*/*.vcf.gz", checkIfExists: !isGCS)
                .filter { vcf -> 
                    vcf.name.endsWith('.vcf.gz') && 
                    !vcf.name.contains('.g.vcf.gz') && 
                    !vcf.name.endsWith('.tbi') 
                }
                .map { vcf -> 
                    def sample = vcf.parent.name
                    tuple(sample, vcf) 
                }
            
            // Collect HaplotypeCaller VCFs
            hc_vcf_ch = Channel
                .fromPath("${variant_calling_outdir}/variant_calling/haplotypecaller/*/*.vcf.gz", checkIfExists: !isGCS)
                .filter { vcf -> 
                    vcf.name.endsWith('.vcf.gz') && 
                    !vcf.name.contains('.g.vcf.gz') && 
                    !vcf.name.endsWith('.tbi') 
                }
                .map { vcf -> 
                    def sample = vcf.parent.name
                    tuple(sample, vcf) 
                }
            
            // Create reference channels
            ref_fasta_ch = Channel.value(file(params.ref_fasta))
            ref_fai_ch = Channel.value(file(params.ref_fasta + ".fai"))
            
            // Run consensus calling
            CONSENSUS_CALLING(dv_vcf_ch, hc_vcf_ch, ref_fasta_ch, ref_fai_ch)
            vcf_ch = CONSENSUS_CALLING.out.consensus_vcf
            
        } else {
            // Use single caller VCF (default to DeepVariant or configurable)
            vcf_ch = Channel
                .fromPath("${variant_calling_outdir}/variant_calling/*/*/*.vcf.gz", checkIfExists: !isGCS)
                .filter { vcf -> 
                    vcf.name.endsWith('.vcf.gz') && 
                    !vcf.name.contains('.g.vcf.gz') && 
                    !vcf.name.endsWith('.tbi') 
                }
                .map { vcf -> 
                    def sample = vcf.parent.name
                    tuple(sample, vcf) 
                }
        }

        // Collect BAMs with BAI
        bam_ch = Channel
            .fromPath("${variant_calling_outdir}/preprocessing/mapped/*/*.sorted.bam", checkIfExists: !isGCS)
            .map { bam -> 
                def sample = bam.parent.name
                def bamPath = bam.toString()
                def baiPath = "${bamPath}.bai"
                
                def bai
                if (isGCS) {
                    bai = file(baiPath, checkIfExists: false)
                } else {
                    bai = file(baiPath)
                    if (!bai.exists()) {
                        bai = file("${bam.parent}/${bam.baseName}.bai")
                        if (!bai.exists()) {
                            error "❌ BAM index not found for ${bam}"
                        }
                    }
                }
                
                tuple(sample, bam, bai) 
            }
        
        // Debug output
        vcf_ch.view { s, v -> "📄 VCF -> ${s} :: ${v.name}" }
        bam_ch.view { s, a, i -> "🧬 BAM -> ${s} :: ${a.name}" }
        
        // Safety checks
        vcf_ch
            .count()
            .subscribe { count ->
                if (count == 0) {
                    error "❌ No VCFs found"
                }
                log.info "✓ Found ${count} VCF file(s)"
            }
        
        bam_ch
            .count()
            .subscribe { count ->
                if (count == 0) {
                    error "❌ No BAMs found in ${variant_calling_outdir}/preprocessing/mapped/*/*.sorted.bam"
                }
                log.info "✓ Found ${count} BAM file(s)"
            }

        // Run post-processing
        POST_SAREK(vcf_ch, bam_ch, bed_ch)
}

workflow RUN_FROM_POST_SAMPLESHEET {
    take:
        post_samplesheet
        bed_ch

    main:
        log.info """
        ╔════════════════════════════════════════════════════════════╗
        ║  Using custom post-samplesheet                             ║
        ║  File: ${post_samplesheet}
        ║  Note: Consensus calling is skipped when using             ║
        ║        post-samplesheet (single VCF per sample expected)   ║
        ╚════════════════════════════════════════════════════════════╝
        """.stripIndent()

        // Parse samplesheet
        Channel
            .fromPath(post_samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def v = file(row.vcf)
                def b = file(row.bam)
                def bi = file(row.bai ?: "${b}.bai")
                
                // Only validate for non-GCS paths
                def isGCS = row.vcf.startsWith('gs://')
                if (!isGCS) {
                    if (!v.exists()) error "❌ VCF not found: ${v}"
                    if (!b.exists()) error "❌ BAM not found: ${b}"
                    if (!bi.exists()) {
                        bi = file("${b.parent}/${b.baseName}.bai")
                        if (!bi.exists()) error "❌ BAI not found for ${b}"
                    }
                }
                
                tuple(row.sample, v, b, bi)
            }
            .multiMap { sample, vcf, bam, bai ->
                vcf: tuple(sample, vcf)
                bam: tuple(sample, bam, bai)
            }
            .set { result }

        vcf_ch = result.vcf
        bam_ch = result.bam
        
        // Debug output
        vcf_ch.view { s, v -> "📄 VCF -> ${s} :: ${v.name}" }
        bam_ch.view { s, a, i -> "🧬 BAM -> ${s} :: ${a.name}" }

        // Run post-processing (no consensus for post-samplesheet)
        POST_SAREK(vcf_ch, bam_ch, bed_ch)
}

workflow RUN_FULL_VARIANT_CALLING {
    take:
        bed_ch

    main:
        log.info """
        ╔════════════════════════════════════════════════════════════╗
        ║  Running full variant calling pipeline (Sarek)             ║
        ║  Consensus calling: ${params.create_consensus ? 'ENABLED' : 'DISABLED'}
        ╚════════════════════════════════════════════════════════════╝
        """.stripIndent()

        PIPELINE_INITIALISATION(
            params.version,
            params.validate_params,
            args,
            params.outdir,
            params.input,
            params.help,
            params.help_full,
            params.show_hidden,
        )
        // Build sample -> assay map from the validated sarek samplesheet channel

        // PIPELINE_INITIALISATION.out.samplesheet is a channel emitting rows (maps/objects)
        // Build sample -> assay map from PIPELINE_INITIALISATION.out.samplesheet
        def assayMap = [:]

        PIPELINE_INITIALISATION.out.samplesheet
            .map { row ->
                def meta = row[0]
                tuple(meta.sample.toString(), (meta.assay ?: 'NA').toString())
            }
            .toList()
            .subscribe { pairs ->
                assayMap = pairs.collectEntries { s, a -> [(s): a] }
                log.info "✓ Loaded assay metadata for ${assayMap.size()} sample(s)"
                assayMap.each { k, v -> log.info "  ${k} -> ${v}" }
            }
        NFCORE_SAREK(PIPELINE_INITIALISATION.out.samplesheet)

        PIPELINE_COMPLETION(
            params.email,
            params.email_on_fail,
            params.plaintext_email,
            params.outdir,
            params.monochrome_logs,
            params.hook_url,
            NFCORE_SAREK.out.multiqc_report
        )

        // Collect variant calling outputs
        COLLECT_VARIANT_CALLING_OUTPUTS(
            NFCORE_SAREK.out.multiqc_report,
            params.outdir
        )

        // Create consensus VCF if enabled
        if (params.create_consensus) {
            ref_fasta_ch = Channel.value(file(params.ref_fasta))
            ref_fai_ch = Channel.value(file(params.ref_fasta + ".fai"))
            
            CONSENSUS_CALLING(
                COLLECT_VARIANT_CALLING_OUTPUTS.out.dv_vcf,
                COLLECT_VARIANT_CALLING_OUTPUTS.out.hc_vcf,
                ref_fasta_ch,
                ref_fai_ch
            )
            
            final_vcf_ch = CONSENSUS_CALLING.out.consensus_vcf
        } else {
            // Use DeepVariant VCFs by default (or could use HC, make configurable)
            final_vcf_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.dv_vcf
        }

        // Run post-processing
        //POST_SAREK(
        //    final_vcf_ch, 
        //    COLLECT_VARIANT_CALLING_OUTPUTS.out.bam, 
        //    bed_ch
        //)
        // Attach meta (sample + assay) to channels for downstream reporting
        def vcf_with_meta_ch = final_vcf_ch.map { sample, vcf ->
            def meta = [ sample: sample, assay: assayMap.get(sample, 'NA') ]
            tuple(meta, vcf)
        }

        def bam_with_meta_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.bam.map { sample, bam, bai ->
            def meta = [ sample: sample, assay: assayMap.get(sample, 'NA') ]
            tuple(meta, bam, bai)
        }


        // Run post-processing (meta-aware)
        POST_SAREK(vcf_with_meta_ch, bam_with_meta_ch, bed_ch)
}

// ═══════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════

workflow {
    
    // Validate and load BED file
    def bedFile = validateBedFile()
    bed_ch = Channel.value(bedFile)
    RUN_FULL_VARIANT_CALLING(bed_ch)
    // Route to appropriate sub-workflow
    //if (params.variant_calling_outdir) {
    //    RUN_FROM_VARIANT_CALLING_OUTDIR(params.variant_calling_outdir, bed_ch)
    //} 
    //else if (params.post_samplesheet) {
    //    RUN_FROM_POST_SAMPLESHEET(params.post_samplesheet, bed_ch)
    //} 
}

// ═══════════════════════════════════════════════════════════════════════════
// WORKFLOW COMPLETION
// ═══════════════════════════════════════════════════════════════════════════

workflow.onComplete {
    log.info """
    ╔════════════════════════════════════════════════════════════╗
    ║  Pipeline completed!                                       ║
    ║  Status: ${workflow.success ? '✓ SUCCESS' : '✗ FAILED'}
    ║  Duration: ${workflow.duration}
    ║  Results: ${params.outdir}
    ╚════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

workflow.onError {
    log.error """
    ╔════════════════════════════════════════════════════════════╗
    ║  ✗ Pipeline failed                                         ║
    ║  Error: ${workflow.errorMessage}
    ╚════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
