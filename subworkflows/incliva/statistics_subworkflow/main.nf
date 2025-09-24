// ===============================
// SUBWORKFLOW: STATISTICS_SUBWORKFLOW
// ===============================

nextflow.enable.dsl=2

include { GLOBAL_MAPPING_STATS } from '../../../modules/incliva/statistics_modules/'

workflow STATISTICS_SUBWORKFLOW {
    take:
    cram_file
    cram_stats
    mosdepth_per_base
    mosdepth_threshold
    fastp_json
    bcftools_stats
    fasta
    fasta_fai
    capture_bed

    main:
    println "[STATISTICS_SUBWORKFLOW] ✅ INICIADO - Procesando estadísticas (modo tumor-only, join version)"

    // Filtrar bcftools → solo mutect2
    ch_bcftools_filtered = bcftools_stats.filter { meta, path -> meta.variantcaller == "mutect2" }

    // Normalizar meta
    ch_cram_file   = cram_file.map        { meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_cram_stats  = cram_stats.map       { meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_per_base    = mosdepth_per_base.map{ meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_thresholds  = mosdepth_threshold.map{ meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_fastp_json  = fastp_json.map       { meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_bcftools    = ch_bcftools_filtered.map{ meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }

    // Join
    ch_joined = ch_cram_file
        .join(ch_cram_stats)
        .join(ch_per_base)
        .join(ch_thresholds)
        .join(ch_fastp_json)
        .join(ch_bcftools)

    // Fasta + Fai
    ch_fasta_pair = fasta.map { meta, f -> f }
                         .combine(fasta_fai.map { meta, f -> f })

    // Capture bed
    ch_bed = capture_bed.map { bed -> bed }

    // Debug
    ch_joined.view { "DEBUG joined → $it" }
    ch_fasta_pair.view { "DEBUG fasta_pair → $it" }
    ch_bed.view { "DEBUG bed → $it" }

    // Llamada al process
    sample_csv = GLOBAL_MAPPING_STATS(ch_joined, ch_fasta_pair, ch_bed)

    // Consolidar resultados en un CSV global
        global_csv = sample_csv.collectFile(
        name: 'global_mapping_stats.csv',
        keepHeader: true,
        storeDir: "${params.outdir}/statistics"
    )


    // Emitir
    emit:
    global_csv
}
