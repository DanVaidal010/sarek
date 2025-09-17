process GLOBAL_MAPPING_STATS {
  publishDir "${params.outdir}/statistics", mode: 'copy'
  label 'global_stats'
  container "/home/danvi/incliva_modulo_estadisticas/singularity/sarek-stats-env.sif"

  input:
  tuple val(meta), path(cram_file)
  tuple val(meta), path(cram_stats)
  tuple val(meta), path(perbase_bed)
  tuple val(meta), path(thresholds_bed)
  tuple val(meta), path(fastp_json)
  tuple val(meta), path(fastq_r1)
  tuple val(meta), path(bcftools_stats)
  tuple val(meta), path(reference_fasta)
  tuple val(meta), path(fasta_fai)
  path capture_bed

  output:
  path "${meta.sample}_mapping_stats.csv"

  script:
  """
  set -euo pipefail


  # 1) Placeholder: resto de métricas de otros scripts
  sample="NA"
  rawdata="NA"
  lowq="NA"
  mapped="NA"
  paired="NA"
  dup="NA"
  mean_insert="NA"
  insert_sd="NA"
  ontarget="NA"
  kitspec="NA"
  cov1q="NA"
  cov3q="NA"
  covmed="NA"
  nt100="NA"
  nt250="NA"
  nt500="NA"
  numvars="NA"

  # 2) Escribir CSV final
  {
    echo "Sample,Rawdata,%LowQReads,%MappedReads,%Paired_reads,%DuplicateReads,Mean_insert,Insert_SD,%OnTargetNoDupReads,Kit_specificity,Cov_1stQ,Cov_3rdQ,Cov_Median,Nt_100x,Nt_250x,Nt_500x,NumVars"
    echo "\$sample,\$rawdata,\$lowq,\$mapped,\$paired,\$dup,\$mean_insert,\$insert_sd,\$ontarget,\$kitspec,\$cov1q,\$cov3q,\$covmed,\$nt100,\$nt250,\$nt500,\$numvars"
  } > ${meta.sample}_mapping_stats.csv
  """
}
