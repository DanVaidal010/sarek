process GLOBAL_MAPPING_STATS {
  tag { meta.patient }
  publishDir "${params.outdir}/${meta.patient}",
             mode: 'copy',
             overwrite: true
  container "${params.container_global_stats}"

  input:
    // Canal con todos los datos de muestra
    tuple val(meta), path(cram), path(stats), path(perbase), path(thr), path(fastp), path(fq1), path(bcftools)
    // Canal con fasta + fai
    tuple path(fasta), path(fai)
    // Canal con capture bed
    path bed

  output:
    path "${meta.patient}_mapping_stats.csv"

  script:
  """
  echo "Sample,CRAM,CRAM_stats,PerBase,Thresholds,FastpJSON,FastqReads,Bcftools,Reference,FastaFAI,CaptureBed" > ${meta.patient}_mapping_stats.csv

  echo "${meta.patient},\\
  ${cram},\\
  ${stats},\\
  ${perbase},\\
  ${thr},\\
  ${fastp},\\
  ${fq1},\\
  ${bcftools},\\
  ${fasta},\\
  ${fai},\\
  ${bed}" >> ${meta.patient}_mapping_stats.csv
  """
}
