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

    // Normalizar meta → usamos una clave única
    ch_cram_file   = cram_file.map        { meta, path -> tuple("${meta.patient}_${meta.sample}", path) }
    ch_cram_stats  = cram_stats.map       { meta, path -> tuple("${meta.patient}_${meta.sample}", path) }
    ch_per_base    = mosdepth_per_base.map{ meta, path -> tuple("${meta.patient}_${meta.sample}", path) }
    ch_thresholds  = mosdepth_threshold.map{ meta, path -> tuple("${meta.patient}_${meta.sample}", path) }
    ch_fastp_json  = fastp_json.map       { meta, path -> tuple("${meta.patient}_${meta.sample}", path) }
    ch_bcftools    = bcftools_stats.filter { meta, path -> meta.variantcaller == "mutect2" }
                                    .map   { meta, path -> tuple("${meta.patient}_${meta.sample}", path) }

    // Debug individuales
    ch_cram_file.view   { "CRAM_FILE KEY → $it" }
    ch_cram_stats.view  { "CRAM_STATS KEY → $it" }
    ch_per_base.view    { "PER_BASE KEY → $it" }
    ch_thresholds.view  { "THRESHOLDS KEY → $it" }
    ch_fastp_json.view  { "FASTP_JSON KEY → $it" }
    ch_bcftools.view    { "BCFTOOLS KEY → $it" }

    // Join por clave única (patient_sample)
    ch_joined = ch_cram_file
        .join(ch_cram_stats)
        .join(ch_per_base)
        .join(ch_thresholds)
        .join(ch_fastp_json)
        .join(ch_bcftools)

    ch_joined.view { "DEBUG joined → $it" }

    // FASTA, FAI y BED como valores constantes
    ch_fasta_val = fasta.map { meta, f -> f }
    ch_fai_val   = fasta_fai.map { meta, f -> f }
    ch_bed_val   = capture_bed


    // Llamada al process
    sample_csv = GLOBAL_MAPPING_STATS(ch_joined, ch_fasta_val, ch_fai_val, ch_bed_val)

    // Consolidar resultados globales
    global_csv = sample_csv.collectFile(
        name: 'global_mapping_stats.csv',
        keepHeader: true,
        storeDir: "${params.outdir}/statistics"
    )

    // Placeholders
    global_amplicon_stats = Channel.empty()
    global_gene_stats     = Channel.empty()

    emit:
    global_csv
    global_amplicon_stats
    global_gene_stats
}
