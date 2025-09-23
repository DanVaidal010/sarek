#!/usr/bin/env bash
set -euo pipefail

# Uso:
#   calculate_kit_specificity_picard.sh <cram> <output_dir|-> <sample_id> <reference_fasta> <capture_bed>

cram_path="$1"
output_dir="$2"
sample_id="$3"
reference_fasta="${4:?missing ref fasta}"
capture_bed="${5:?missing capture bed}"

prefix="${output_dir}/${sample_id}"
bam="${prefix}.bam"
bam_target="${prefix}.ontarget.bam"
metrics_all="${prefix}_mapped_alignment_metrics.txt"
metrics_target="${prefix}_ontarget_alignment_metrics.txt"
result_file="${prefix}.kit_specificity_result.txt"

cleanup() {
  [[ "$output_dir" == "-" ]] && {
    rm -f "$bam" "$bam.bai" "$bam_target" "$bam_target.bai" "$metrics_all" "$metrics_target" 2>/dev/null || true
  }
}
trap cleanup EXIT

if [[ "$output_dir" == "-" ]]; then
  prefix="${sample_id}.tmp"
  bam="${prefix}.bam"
  bam_target="${prefix}.ontarget.bam"
  metrics_all="${prefix}_mapped_alignment_metrics.txt"
  metrics_target="${prefix}_ontarget_alignment_metrics.txt"
  result_file="${prefix}.kit_specificity_result.txt"
fi

# Convertir CRAM a BAM
samtools view -b -T "$reference_fasta" "$cram_path" > "$bam"
samtools index "$bam"

# Filtrar BAM para solo on-target
samtools view -b -L "$capture_bed" "$bam" > "$bam_target"
samtools index "$bam_target"

# Picard para métricas
picard CollectAlignmentSummaryMetrics \
    I="$bam" \
    O="$metrics_all" \
    R="$reference_fasta" \
    VALIDATION_STRINGENCY=LENIENT

picard CollectAlignmentSummaryMetrics \
    I="$bam_target" \
    O="$metrics_target" \
    R="$reference_fasta" \
    VALIDATION_STRINGENCY=LENIENT

# Calcular Kit Specificity
mapped_reads=$(grep -m1 "^PAIR" "$metrics_all" | awk '{print $6}' || echo "0")
ontarget_reads=$(grep -m1 "^PAIR" "$metrics_target" | awk '{print $6}' || echo "0")

if [[ -z "$mapped_reads" || -z "$ontarget_reads" || "$mapped_reads" -eq 0 ]]; then
  specificity="0.00"
else
  specificity=$(awk "BEGIN { printf \"%.2f\", ($ontarget_reads / $mapped_reads) * 100 }")
fi

# Output
if [[ "$output_dir" == "-" ]]; then
  echo "$specificity"
else
  echo -e "Kit_specificity\t$specificity" > "$result_file"
  echo "✅ Resultado guardado en $result_file"
fi
