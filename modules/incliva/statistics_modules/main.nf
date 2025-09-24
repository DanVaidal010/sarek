process GLOBAL_MAPPING_STATS {
    tag "${key}"

    publishDir "${params.outdir}/${key}",
        mode: 'copy',
        overwrite: true

    container "${params.container_global_stats}"

    input:
    tuple val(key), path(cram), path(stats), path(perbase), path(thr), path(fastp), path(bcftools)
    val fasta
    val fai
    path bed

    output:
    path "${key}_mapping_stats.csv"

    script:
    """
    echo "▶️ [GLOBAL_MAPPING_STATS] Procesando muestra ${key}" >&2

    # === Definir archivos intermedios ===
    basic_output="${key}.fastq_summary.txt"
    cuartil_output="${key}.coverage_quartiles.txt"
    mean_output="${key}.coverage_thresholds.txt"
    insert_output="${key}.insert_metrics.txt"
    numvars_output="${key}.variant_count.txt"
    ontarget_output="${key}.ontarget_result.txt"
    kit_output="${key}.kit_specificity_result.txt"

    # 1. Estadísticas básicas (fastp + cram.stats)
    summarize_fastq_cram_stats.sh ${fastp} ${stats} ${key} \$basic_output

    # 2. Cuartiles de cobertura
    extract_coverage_quartiles.py -p ${perbase} -c ${bed} -o \$cuartil_output

    # 3. Cobertura ≥100x/250x/500x
    coverage_depth_thresholds_100_250_500.py -i ${thr} -o \$mean_output

    # 4. Inserto medio y desviación estándar
    calculate_insert_metrics.py -j ${fastp} -o \$insert_output

    # 5. Conteo de variantes
    count_filtered_variants.sh ${bcftools} \$numvars_output

    # 6. %OnTargetNoDupReads con Picard
    calculate_ontarget_picard.sh ${cram} ${fastp} ${fasta} ${key} .

    # 7. Kit Specificity con Picard
    calculate_kit_specificity_picard.sh ${cram} ${fasta} ${bed} ${key} .

    # 8. Combinar resultados en CSV final
    combine_mapping_statistics.py \\
        --sample ${key} \\
        --quartile \$cuartil_output \\
        --mean \$mean_output \\
        --insert \$insert_output \\
        --numvars \$numvars_output \\
        --basic \$basic_output \\
        --ontarget \$ontarget_output \\
        --kit \$kit_output \\
        --output ${key}_mapping_stats.csv
    """
}
