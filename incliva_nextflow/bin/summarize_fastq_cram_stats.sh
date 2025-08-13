#!/bin/bash
set -euo pipefail

# Uso:
# ./extract_basic_stats.sh fastp.json cram.stats sample_name output_file

json_file="$1"      # Archivo JSON de fastp
stats_file="$2"     # Archivo .cram.stats de samtools
sample="$3"         # Nombre de la muestra
output_file="$4"    # Ruta del archivo de salida

# Extraer datos desde JSON usando jq
passed=$(jq '.filtering_result.passed_filter_reads' "$json_file")
lowq=$(jq '.filtering_result.low_quality_reads' "$json_file")
dup_rate=$(jq '.duplication.rate' "$json_file")

# Extraer estadísticas básicas desde cram.stats
total=$(grep "raw total sequences:" "$stats_file" | awk -F':' '{print $2}' | awk '{print $1}')
mapped=$(grep "reads mapped:" "$stats_file" | awk -F':' '{print $2}' | awk '{print $1}')
paired=$(grep "reads paired:" "$stats_file" | awk -F':' '{print $2}' | awk '{print $1}')
properly_paired=$(grep "reads properly paired:" "$stats_file" | awk -F':' '{print $2}' | awk '{print $1}')

# Calcular % de lecturas de baja calidad
lowq_pct=$(awk "BEGIN {printf \"%.2f\", 100 * $lowq / ($passed + $lowq)}")

# Forzar % de duplicados a 0.00 (en caso de amplicones)
dup_pct=0.00

# Calcular % mapeadas y % correctamente emparejadas
mapped_pct=$(awk "BEGIN {printf \"%.2f\", 100 * $mapped / $total}")
paired_pct=$(awk "BEGIN {printf \"%.2f\", 100 * $properly_paired / $paired}")

# Escribir resultados al archivo
{
  echo -e "Sample\t$sample"
  echo -e "Rawdata\t$passed"
  echo -e "%LowQReads\t$lowq_pct"
  echo -e "%MappedReads\t$mapped_pct"
  echo -e "%Paired_reads\t$paired_pct"
  echo -e "%DuplicateReads\t$dup_pct"
} > "$output_file"

echo "✅ Guardado en $output_file"
