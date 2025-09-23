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
    fastq_reads
    bcftools_stats
    fasta
    fasta_fai
    capture_bed

    main:
    //
    // 1. Filtrar bcftools → solo mutect2
    //
    ch_bcftools_filtered = bcftools_stats.filter { meta, path ->
        meta.variantcaller == "mutect2"
    }

    //
    // 2. Normalizar todos los meta (solo patient y sample)
    //
    ch_cram_file   = cram_file.map        { meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_cram_stats  = cram_stats.map       { meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_per_base    = mosdepth_per_base.map{ meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_thresholds  = mosdepth_threshold.map{ meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_fastp_json  = fastp_json.map       { meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_fastq_reads = fastq_reads.map      { meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }
    ch_bcftools    = ch_bcftools_filtered.map{ meta, path -> tuple([patient: meta.patient, sample: meta.sample], path) }

    //
    // 3. Unir canales por sample
    //
    ch_joined = ch_cram_file
        .join(ch_cram_stats)
        .join(ch_per_base)
        .join(ch_thresholds)
        .join(ch_fastp_json)
        .join(ch_fastq_reads)
        .join(ch_bcftools)

    //
    // 4. Fasta y Fai sin meta (nos quedamos solo con los paths)
    //
    ch_fasta_pair = fasta.map { meta, f -> f }
                         .combine(fasta_fai.map { meta, f -> f })

    //
    // 5. Capture bed es un path único → lo pasamos limpio
    //
    ch_bed = capture_bed.map { bed -> bed }

    //
    // 6. Debug para ver qué llega
    //
    ch_joined.view { "DEBUG joined → $it" }
    ch_fasta_pair.view { "DEBUG fasta_pair → $it" }
    ch_bed.view { "DEBUG bed → $it" }

    //
    // 7. Llamada al process
    //
    GLOBAL_MAPPING_STATS(ch_joined, ch_fasta_pair, ch_bed)

    emit:
    global_csv = GLOBAL_MAPPING_STATS.out
}
