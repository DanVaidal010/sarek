nextflow.enable.dsl=2

include { GLOBAL_MAPPING_STATS } from '../../../modules/incliva/statistics_modules/'

workflow STATISTICS_SUBWORKFLOW {

  take:
    input_sample
    cram_file
    cram_stats
    mosdepth_per_base
    mosdepth_threshold
    fastp_json
    fastq_reads
    bcftools_stats
    fasta
    fasta_fai
    capture_bed

  main:
    println "[STATISTICS_SUBWORKFLOW] ✅ INICIADO - Procesando estadísticas (modo tumor-only)"

    if (params.debug_statistics) {
        input_sample.view { meta_and_files ->
            def m = meta_and_files[0]
            "[DEBUG] (subworkflow) Sample=${m.sample}, Status=${m.status}"
        }

        // Debug de canales principales
        cram_file.view      { "[DEBUG] cram_file → ${it}" }
        cram_stats.view     { "[DEBUG] cram_stats → ${it}" }
        fastp_json.view     { "[DEBUG] fastp_json → ${it}" }
        bcftools_stats.view { "[DEBUG] bcftools_stats → ${it}" }
        fasta.view          { "[DEBUG] fasta → ${it}" }
        fasta_fai.view      { "[DEBUG] fasta_fai → ${it}" }
        capture_bed.view    { "[DEBUG] capture_bed → ${it}" }
    }

    // ⚠️ Solo normalizamos fastq_reads para quedarnos con R1
    ch_fastq_r1 = fastq_reads.map { meta, files ->
        def r1 = (files instanceof List && files) ? files[0] : files
        tuple(meta, r1)
    }

    // Llamar al módulo GLOBAL_MAPPING_STATS
    sample_csv = GLOBAL_MAPPING_STATS(
        cram_file,
        cram_stats,
        mosdepth_per_base,
        mosdepth_threshold,
        fastp_json,
        ch_fastq_r1,
        bcftools_stats,
        fasta,
        fasta_fai,
        capture_bed
    )

    // Consolidar resultados globales
    global_csv = sample_csv.collectFile(
        name: 'global_mapping_stats.csv',
        keepHeader: true,
        storeDir: "${params.outdir}/statistics"
    )

    // Placeholders para futuras métricas
    global_amplicon_stats = Channel.empty()
    global_gene_stats     = Channel.empty()

  emit:
    global_csv
    global_amplicon_stats
    global_gene_stats
}
