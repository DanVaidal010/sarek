#!/bin/bash

# Uso:
# ./extract_basic_stats.sh fastp.json cram.stats sample_name output_file

json_file="$1"
stats_file="$2"
sample="$3"
output_file="$4"

# Extraer valores de fastp.json con jq
passed=$(jq '.filtering_result.passed_filter_reads' "$json_file")
lowq=$(jq '.filtering_result.low_quality_reads' "$json_file")
dup_rate=$(jq '.duplication.rate' "$json_file")

# Extraer valores del archivo .cram.stats con grep y awk
total=$(grep "raw total sequences:" "$stats_file" | awk -F':' '{print $2}' | awk '{print $1}')
mapped=$(grep "reads mapped:" "$stats_file" | awk -F':' '{print $2}' | awk '{print $1}')
paired=$(grep "reads paired:" "$stats_file" | awk -F':' '{print $2}' | awk '{print $1}')
properly_paired=$(grep "reads properly paired:" "$stats_file" | awk -F':' '{print $2}' | awk '{print $1}')

# Cálculos
lowq_pct=$(awk "BEGIN {printf \"%.2f\", 100 * $lowq / ($passed + $lowq)}")

# 🔧 Forzar duplicados a 0 si se trabaja con amplicones
dup_pct=0.00
#dup_pct=$(awk "BEGIN {printf \"%.2f\", 100 * $dup_rate}")

mapped_pct=$(awk "BEGIN {printf \"%.2f\", 100 * $mapped / $total}")
paired_pct=$(awk "BEGIN {printf \"%.2f\", 100 * $properly_paired / $paired}")

# Escribir salida (sin incluir la línea "Sample\t<sample>")
{
  echo -e "Sample\t$sample" > "$output_file"
  echo -e "Rawdata\t$passed"
  echo -e "%LowQReads\t$lowq_pct"
  echo -e "%MappedReads\t$mapped_pct"
  echo -e "%Paired_reads\t$paired_pct"
  echo -e "%DuplicateReads\t$dup_pct"
} > "$output_file"


echo "✅ Guardado en $output_file"

#Poner el duplicates =0