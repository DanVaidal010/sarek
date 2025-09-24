process GLOBAL_MAPPING_STATS {
    tag "${meta.patient}_${meta.sample}"

    publishDir "${params.outdir}/${meta.sample}",
        mode: 'copy',
        overwrite: true

    container "${params.container_global_stats}"

    input:
    tuple val(meta), path(cram), path(stats), path(perbase), path(thr), path(fastp), path(bcftools)
    tuple path(fasta), path(fai)
    path bed

    output:
    path "${meta.sample}_mapping_stats.csv"

    script:
    """
    # Cabecera
    echo "Sample,CRAM,CRAM_stats,PerBase,Thresholds,FastpJSON,Bcftools,Reference,FastaFAI,CaptureBed" > ${meta.sample}_mapping_stats.csv

    # Fila con nombres de archivo
    echo "${meta.sample},\\
    \$(basename ${cram}),\\
    \$(basename ${stats}),\\
    \$(basename ${perbase}),\\
    \$(basename ${thr}),\\
    \$(basename ${fastp}),\\
    \$(basename ${bcftools}),\\
    \$(basename ${fasta}),\\
    \$(basename ${fai}),\\
    \$(basename ${bed})" >> ${meta.sample}_mapping_stats.csv
    """
}
