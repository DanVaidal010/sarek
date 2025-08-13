nextflow.enable.dsl = 2

// Usa params.modules_dir en vez de rutas absolutas
include { GLOBAL_MAPPING_STATS } from "${params.modules_dir}/statistics_modules.nf"
// Cuando integremos amplicón/gen, descomenta estos includes y los bloques del main
// include { STATS_BY_AMPLICON               } from "${params.modules_dir}/statistics_modules.nf"
// include { CREATE_GLOBAL_STATS_BY_AMPLICON } from "${params.modules_dir}/statistics_modules.nf"
// include { STATS_BY_GENE                   } from "${params.modules_dir}/statistics_modules.nf"
// include { CREATE_GLOBAL_STATS_BY_GENE     } from "${params.modules_dir}/statistics_modules.nf"

workflow STATISTICS_SUBWORKFLOW {

  // ----------- ENTRADAS (mismo orden que en workflows/sarek/main.nf) -----------
  take:
    ch_cram_pair_to_cross          // (patient, meta, cram, crai)
    ch_intervals_and_num_intervals // (BED, n) — no usado por GLOBAL_MAPPING_STATS
    ch_capture_amplicon            // value channel: path(panel.bed)
    val outdir                     // no usado (el process usa params.outdir)
    ch_fasta                       // value channel: path(fasta)
    ch_fasta_fai                   // value channel: path(fasta.fai)
    ch_samtools_stats              // tuple val(meta), path(samtools_stats.txt)
    ch_fastp_json                  // tuple val(meta), path(fastp.json)
    ch_bcftools_stats              // tuple val(meta), path(bcftools_stats.txt)
    ch_mosdepth_threshold          // tuple val(meta), path(thresholds.bed.gz)
    ch_mosdepth_per_base           // tuple val(meta), path(per-base.bed.gz)

  // ---------------- CUERPO ----------------
  main:
    // 1) Adaptar el CRAM: (patient, meta, cram, crai) -> (meta, cram)
    ch_cram = ch_cram_pair_to_cross.map { patient, meta, cram, crai ->
      tuple(meta, cram)
    }

    // 2) Resto de canales Sarek-style tal cual
    ch_sam   = ch_samtools_stats
    ch_fp    = ch_fastp_json
    ch_bcf   = ch_bcftools_stats
    ch_thr   = ch_mosdepth_threshold
    ch_pbase = ch_mosdepth_per_base

    // 3) Broadcast de FASTA/FAI/BED a cada muestra usando los value channels de entrada
    ch_ref = ch_cram
      .map { meta, _ -> meta }           // nos quedamos con 'meta'
      .combine(ch_fasta)                 // (meta, fasta)
      .combine(ch_fasta_fai)             // (meta, fasta, fai)
      .combine(ch_capture_amplicon)      // (meta, fasta, fai, bed)
      .map { meta, fasta, fai, bed -> tuple(meta, fasta, fai, bed) }

    // 4) Invocar el process y recoger su canal de salida con .out
    GLOBAL_MAPPING_STATS(
      ch_cram,      // (meta, cram)
      ch_sam,       // (meta, samtools_stats)
      ch_pbase,     // (meta, per-base.bed.gz)
      ch_thr,       // (meta, thresholds.bed.gz)
      ch_fp,        // (meta, fastp.json)
      ch_bcf,       // (meta, bcftools_stats)
      ch_ref        // (meta, fasta, fai, bed)
    )
    mapping_csvs = GLOBAL_MAPPING_STATS.out

    // 5) CSV global agregado — patrón Sarek
    global_mapping_stats = mapping_csvs.collectFile(
      name: "${params.study_name ?: 'study'}_global_mapping_stats.csv",
      sort: true,
      keepHeader: true,
      storeDir: "${params.outdir}"
    )

    // ----------- ESTADÍSTICAS POR AMPLICÓN / GEN -----------
    // Aún no integradas: emitimos canales vacíos para no romper el .mix() del caller
    global_amplicon_stats = Channel.empty()
    global_gene_stats     = Channel.empty()

  // ----------- SALIDAS DEL SUBWORKFLOW -----------
  emit:
    global_mapping_stats
    global_amplicon_stats
    global_gene_stats
}
