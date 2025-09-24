process GLOBAL_MAPPING_STATS {
    tag "${key}"  // tag único por muestra

    publishDir "${params.outdir}/${key}",
        mode: 'copy',
        overwrite: true

    container "${params.container_global_stats}"

    input:
    tuple val(key), path(cram), path(stats), path(perbase), path(thr), path(fastp), path(bcftools)
    val fasta
    val fai
    val bed

    output:
    path "${key}_mapping_stats.csv"

    script:
    """
    echo "DEBUG → Running GLOBAL_MAPPING_STATS for sample: ${key}" >&2

    echo "Sample,CRAM,CRAM_stats,PerBase,Thresholds,FastpJSON,Bcftools,Reference,FastaFAI,CaptureBed" > ${key}_mapping_stats.csv

    echo "${key},\\
    \$(basename ${cram}),\\
    \$(basename ${stats}),\\
    \$(basename ${perbase}),\\
    \$(basename ${thr}),\\
    \$(basename ${fastp}),\\
    \$(basename ${bcftools}),\\
    \$(basename ${fasta}),\\
    \$(basename ${fai}),\\
    \$(basename ${bed})" >> ${key}_mapping_stats.csv
    """
}
