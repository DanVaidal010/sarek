process GLOBAL_MAPPING_STATS {                                  // Declara un proceso Nextflow (unidad de ejecución)

  tag { meta.sample }                                           // Etiqueta de tarea: se muestra en logs/monitor por muestra
  publishDir "${params.outdir}/${meta.sample}",                // Directorio donde se publica el resultado final
            mode: 'copy',                                      // Copia el/los archivo(s) de salida
            overwrite: true                                    // Sobrescribe si ya existían
  container "${params.container_global_stats}"                 // Contenedor con samtools, jq, python3 y picard (ajusta en params)

  input:                                                       // Definición de entradas (canales Sarek-style)
  tuple val(meta), path(cram_file)                             // 1) CRAM recalibrado de la muestra
  tuple val(meta), path(cram_stats)                            // 2) Stats de samtools (samtools stats) del CRAM
  tuple val(meta), path(perbase_bed)                           // 3) Mosdepth per-base (*.per-base.bed.gz) para cuantiles
  tuple val(meta), path(thresholds_bed)                        // 4) Mosdepth thresholds (*.thresholds.bed.gz) 100/250/500x
  tuple val(meta), path(fastp_json)                            // 5) JSON de fastp (lowQ, insert histogram, totales, etc.)
  tuple val(meta), path(fastq_r1, optional: true)              // 6) [Opcional] FASTQ R1 (si existe): para usar como denominador “raw”
  tuple val(meta), path(vcf_file,  optional: true)             // 7) [Opcional] VCF filtrado por muestra (si no, se usa bcftools_stats)
  tuple val(meta), path(bcftools_stats, optional: true)        // 8) [Opcional] bcftools stats (si está, mejor para num_vars)
  tuple val(meta), path(reference_fasta),                      // 9) FASTA de referencia (para Picard/samtools view)
                 path(reference_fai),                          //    FAI del FASTA
                 path(capture_bed)                             //    Panel de captura (BED) para cuantiles y Kit specificity

  output:                                                      // Salidas del proceso
  path "${meta.sample}_mapping_stats.csv"                      // CSV único por muestra (nombre determinista)

  script:                                                      // Script de la tarea (bash + tus scripts .py)
  """
  set -euo pipefail                                            # Fail rápido en errores/variables no definidas/pipes

  # ------------------------------------------------------------------------------
  # Helpers Picard (funciones bash)                                          #
  # ------------------------------------------------------------------------------

  # calculate_usable_reads_picard:
  # - Convierte CRAM→BAM, corre Picard CollectAlignmentSummaryMetrics
  # - Toma PF_READS_ALIGNED (col 6, fila 'PAIR') como “lecturas alineadas útiles”
  # - Denominador: si hay FASTQ R1 -> líneas/2 ; si no hay -> total_reads de fastp.json
  
  calculate_usable_reads_picard() {
    local cram_path="\$1" ref="\$2" fastq_r1_opt="\$3" fastp_json="\$4"  # Args: CRAM, FASTA, [R1 opcional], fastp.json
    local bam_raw="\${meta.sample}.usable.tmp.bam"                        # BAM temporal
    local metrics="\${meta.sample}.usable.metrics.txt"                    # Métricas Picard temporales

    samtools view -b -T "\$ref" "\$cram_path" > "\$bam_raw"              # CRAM -> BAM usando la referencia
    samtools index "\$bam_raw" >/dev/null                                # Indexa el BAM

    picard CollectAlignmentSummaryMetrics \                               # Ejecuta Picard para métricas globales
      I="\$bam_raw" O="\$metrics" R="\$ref" VALIDATION_STRINGENCY=LENIENT

    local aligned_reads                                                   # PF_READS_ALIGNED (aprox. col 6 en línea 'PAIR')
    aligned_reads=\$(grep -m1 "^PAIR" "\$metrics" | awk '{print \$6}' || echo "0")

    local total_reads="0"                                                 # Denominador (total de lecturas)
    if [[ -n "\${fastq_r1_opt:-}" && -f "\$fastq_r1_opt" ]]; then        # Si tenemos R1, usamos ese (líneas/2)
      total_reads=\$(zcat "\$fastq_r1_opt" | wc -l | awk '{printf "%d", \$1/2}')
    else                                                                  # Si NO hay R1, usamos fastp.json (total_reads)
      total_reads=\$(jq -r '.summary.before_filtering.total_reads // (.filtering_result.passed_filter_reads + .filtering_result.low_quality_reads)' "\$fastp_json" 2>/dev/null || echo "0")
    fi

    local pct="0.00"                                                      # Porcentaje de “usable reads”
    if [[ "\$aligned_reads" =~ ^[0-9]+$ && "\$total_reads" =~ ^[0-9]+$ && "\$total_reads" -gt 0 ]]; then
      pct=\$(awk "BEGIN { printf \\"%.2f\\", ( \$aligned_reads / \$total_reads ) * 100 }")  # aligned / total * 100
    fi

    rm -f "\$bam_raw" "\$bam_raw.bai" "\$metrics"                         # Limpieza de temporales
    echo "\$pct"                                                          # Devuelve el % por stdout
  }

  # calculate_kit_specificity:
  # - CRAM→BAM, filtra on-target con -L capture.bed, corre Picard en ambos
  # - (lecturas mapeadas on-target / mapeadas global) * 100
  calculate_kit_specificity() {
    local cram_path="\$1" ref="\$2" bed="\$3"                             # Args: CRAM, FASTA, BED panel
    local bam="\${meta.sample}.kit.tmp.bam"                               # BAM global temporal
    local bam_target="\${meta.sample}.kit.tmp.ontarget.bam"               # BAM on-target temporal
    local m_all="\${meta.sample}.kit.tmp_mapped_all.txt"                  # Métricas Picard global
    local m_tar="\${meta.sample}.kit.tmp_ontarget_metrics.txt"            # Métricas Picard on-target

    samtools view -b -T "\$ref" "\$cram_path" > "\$bam"                   # CRAM -> BAM
    samtools index "\$bam" >/dev/null                                     # Index global
    samtools view -b -L "\$bed" "\$bam" > "\$bam_target"                  # Filtra on-target por panel
    samtools index "\$bam_target" >/dev/null                              # Index on-target

    picard CollectAlignmentSummaryMetrics I="\$bam"        O="\$m_all" R="\$ref" VALIDATION_STRINGENCY=LENIENT   # Métricas globales
    picard CollectAlignmentSummaryMetrics I="\$bam_target" O="\$m_tar" R="\$ref" VALIDATION_STRINGENCY=LENIENT   # Métricas on-target

    local mapped ontarget                                                 # Extrae lecturas mapeadas
    mapped=\$(  grep -m1 "^PAIR" "\$m_all" | awk '{print \$6}' || echo "0")
    ontarget=\$(grep -m1 "^PAIR" "\$m_tar"  | awk '{print \$6}' || echo "0")

    local spec="0.00"                                                     # % especificidad
    if [[ "\$mapped" =~ ^[0-9]+$ && "\$ontarget" =~ ^[0-9]+$ && "\$mapped" -gt 0 ]]; then
      spec=\$(awk "BEGIN { printf \\"%.2f\\", ( \$ontarget / \$mapped ) * 100 }")           # on-target / mapeadas
    fi

    rm -f "\$bam" "\$bam.bai" "\$bam_target" "\$bam_target.bai" "\$m_all" "\$m_tar"          # Limpieza
    echo "\$spec"                                                         # Devuelve el % por stdout
  }

  # ------------------------------------------------------------------------------
  # 1) Métricas rápidas de samtools & fastp                                 #
  # ------------------------------------------------------------------------------

  total=\$(grep -m1 "raw total sequences:"               "${cram_stats}" | awk -F':' '{print \$2}' | awk '{print \$1}')  # Total lecturas
  mapped=\$(grep -m1 "reads mapped:"                     "${cram_stats}" | awk -F':' '{print \$2}' | awk '{print \$1}')  # Mapeadas
  paired=\$(grep -m1 "reads paired:"                     "${cram_stats}" | awk -F':' '{print \$2}' | awk '{print \$1}')  # Pareadas
  properly_paired=\$(grep -m1 "reads properly paired:"   "${cram_stats}" | awk -F':' '{print \$2}' | awk '{print \$1}')  # Pareadas “properly”

  passed=\$(jq '.filtering_result.passed_filter_reads'    "${fastp_json}")                                            # fastp: pasadas filtro
  lowq=\$(jq '.filtering_result.low_quality_reads'        "${fastp_json}")                                            # fastp: baja calidad

  raw_data="\$total"                                                                                                  # Guardamos total
  lowq_pct=\$(awk "BEGIN {p=\$passed+\$lowq; if(p>0) printf \\"%.2f\\", 100*\$lowq/p; else print \\"0.00\\" }")       # % lowQ
  mapped_pct=\$(awk "BEGIN { if(\$total>0)   printf \\"%.2f\\", 100*\$mapped/\$total;   else print \\"0.00\\" }")     # % mapeadas
  paired_pct=\$(awk "BEGIN { if(\$paired>0)  printf \\"%.2f\\", 100*\$properly_paired/\$paired; else print \\"0.00\\" }") # % properly paired
  dup_pct="0.00"                                                                                                     # Amplicones: sin duplicados

  # ------------------------------------------------------------------------------
  # 2) Cuantiles de cobertura (per-base + panel)                            #
  # ------------------------------------------------------------------------------

  quart_out=\$(python3 extract_coverage_quartiles.py -p "${perbase_bed}" -c "${capture_bed}" -o -)   # Q1/Q3/mediana al stdout
  q1=\$(  printf "%s\\n" "\$quart_out" | awk 'NR==1{print \$2}')                                      # Q1
  q3=\$(  printf "%s\\n" "\$quart_out" | awk 'NR==2{print \$2}')                                      # Q3
  qmed=\$(printf "%s\\n" "\$quart_out" | awk 'NR==3{print \$2}')                                      # Mediana

  # ------------------------------------------------------------------------------
  # 3) Cobertura a umbrales 100x/250x/500x                                  #
  # ------------------------------------------------------------------------------

  thr_out=\$(python3 coverage_depth_thresholds_100_250_500.py -i "${thresholds_bed}" -o -)            # Procesa thresholds.bed.gz
  d100=\$(printf "%s\\n" "\$thr_out" | awk 'NR==1{print \$2}')                                         # Nt_100x
  d250=\$(printf "%s\\n" "\$thr_out" | awk 'NR==2{print \$2}')                                         # Nt_250x
  d500=\$(printf "%s\\n" "\$thr_out" | awk 'NR==3{print \$2}')                                         # Nt_500x

  # ------------------------------------------------------------------------------
  # 4) Inserto medio y desviación estándar                                  #
  # ------------------------------------------------------------------------------

  ins_out=\$(python3 calculate_insert_metrics.py -j "${fastp_json}" -o -)                              # Usa histograma de fastp.json
  ins_mean=\$(printf "%s\\n" "\$ins_out" | awk '/Mean_insert/{print \$2}')                             # insert_mean
  ins_sd=\$(  printf "%s\\n" "\$ins_out" | awk '/Insert_SD/{print \$2}')                               # insert_SD

  # ------------------------------------------------------------------------------
  # 5) %Usable reads y Kit specificity                                     #
  # ------------------------------------------------------------------------------

  usable=\$(calculate_usable_reads_picard  "${cram_file}" "${reference_fasta}" "${fastq_r1:-}" "${fastp_json}")  # %usable_reads
  kitsp=\$( calculate_kit_specificity      "${cram_file}" "${reference_fasta}" "${capture_bed}")                  # Kit_specificity

  # ------------------------------------------------------------------------------
  # 6) Número de variantes (preferencia: bcftools_stats → VCF → NA)        #
  # ------------------------------------------------------------------------------

  if [[ -s "${bcftools_stats:-}" ]]; then                                                              # Si tenemos bcftools_stats
    numvars=\$(awk '/^SN/ && \$3=="records:" {print \$4}' "${bcftools_stats}" | tail -n1)              # Toma “number of records”
    numvars=\${numvars:-NA}                                                                            # Si vacío, pone NA
  elif [[ -s "${vcf_file:-}" ]]; then                                                                  # Si no, intenta VCF
    numvars=\$(grep -vc "^#" "${vcf_file}")                                                            # Cuenta no-cabeceras
  else
    numvars="NA"                                                                                       # Si nada, NA
  fi

  # ------------------------------------------------------------------------------
  # 7) CSV final único por muestra                                         #
  # ------------------------------------------------------------------------------

  {
  echo "Sample,Rawdata,%LowQReads,%MappedReads,%Paired_reads,%DuplicateReads,%OnTargetNoDupReads,Kit_specificity,Cov_1stQ,Cov_3rdQ,Cov_Median,Nt_100x,Nt_250x,Nt_500x,Mean_insert,Insert_SD,NumVars"
  echo "${meta.sample},${raw_data},${lowq_pct},${mapped_pct},${paired_pct},${dup_pct},${usable},${kitsp},${q1},${q3},${qmed},${d100},${d250},${d500},${ins_mean},${ins_sd},${numvars}"
} > "${meta.sample}_mapping_stats.csv"                                                                # Escribe CSV final
  """
}

process STATS_BY_GENE {
  tag { meta.sample }                                                   // etiqueta por muestra
  publishDir "${params.outdir}/${meta.sample}", mode:'copy', overwrite:true
  container "${params.container_global_stats}"                          // contenedor con python3

  input:
  tuple val(meta), path(perbase_bed)                                    // mosdepth *.per-base.bed.gz
  tuple val(meta), path(thresholds_bed)                                 // reservado (no se usa en el script, pero útil si amplías)
  tuple val(meta), path(capture_bed)                                    // BED del panel (chr start end [amplicon_id] [gene])
  tuple val(meta), path(annotated_csv, optional: true)                  // annotated.csv (opcional)

  output:
  tuple val(meta), path("${meta.sample}.gene_stats.tsv"), emit: gene_tsv // TSV por muestra (cabecera Excel)

  /*
   * Notas:
   * - params.coverage_cutoff define el umbral (p.ej. 250x) para "covered".
   * - Si NO pasas annotated.csv, el script pone "0" y "-" en las dos columnas de variantes.
   */
  script:
  """
  set -euo pipefail

  compute_target_coverage.py \\
    --mode gene \\
    --per-base "${perbase_bed}" \\
    --capture-bed "${capture_bed}" \\
    --cutoff "${params.coverage_cutoff ?: 250}" \\
    ${annotated_csv ? "--annotated ${annotated_csv}" : ""} \\
    --sample "${meta.sample}" \\
    --output "${meta.sample}.gene_stats.tsv"
  """
}


process CREATE_GLOBAL_STATS_BY_GENE {
  tag "global_gene"
  publishDir "${params.outdir}", mode:'copy', overwrite:true
  container "${params.container_global_stats}"   // mismo contenedor con python3

  input:
  path(gene_tsv_list)                            // viene de gene_tsv.collect()

  output:
  path("${params.study_name}_stats_by_gene.tsv"), emit: global_gene_stats

  /*
   * Une los TSV por muestra en una tabla global ancha por GENE.
   * Cabecera final: GENE + bloques de 3 columnas por muestra (ver script).
   */
  script:
  """
  set -euo pipefail
  merge_gene_stats.py \\
    --inputs ${gene_tsv_list.join(' ')} \\
    --output "${params.study_name}_stats_by_gene.tsv"
  """
}

