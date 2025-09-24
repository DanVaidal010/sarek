#!/usr/bin/env bash
set -euo pipefail

# === INPUTS ===
cram_path="$1"        # CRAM recalibrado
fastp_json="$2"       # fastp.json
reference_fasta="$3"  # FASTA de referencia
sample_id="$4"        # ID de muestra
output_dir="$5"       # Carpeta de salida

# === OUTPUTS ===
prefix="${output_dir}/${sample_id}"
bam_raw="${prefix}.mapped.bam"
metrics_file="${prefix}_ontarget_alignment_metrics.txt"
result_file="${prefix}.ontarget_result.txt"

echo "🧪 Procesando %OnTargetNoDupReads para muestra: $sample_id"

# 1. CRAM → BAM
samtools view -b -T "$reference_fasta" "$cram_path" > "$bam_raw"
samtools index "$bam_raw"

# 2. Ejecutar Picard
picard CollectAlignmentSummaryMetrics \
    I="$bam_raw" \
    O="$metrics_file" \
    R="$reference_fasta" \
    VALIDATION_STRINGENCY=LENIENT

# 3. Extraer lecturas on-target
ontarget_reads=$(grep -m1 "^PAIR" "$metrics_file" | awk '{print $6}' || echo "0")

# 4. Extraer lecturas crudas de fastp.json
raw_reads=$(jq -r '.summary.before_filtering.total_reads' "$fastp_json")

# 5. Calcular porcentaje
if [[ -z "$ontarget_reads" || -z "$raw_reads" || "$raw_reads" -eq 0 ]]; then
    on_target_pct="0.00"
else
    on_target_pct=$(awk "BEGIN { printf \"%.2f\", ($ontarget_reads / $raw_reads) * 100 }")
fi

# 6. Guardar resultado
echo -e "%OnTargetNoDupReads\t$on_target_pct" > "$result_file"
echo "✅ Resultado guardado en $result_file"

# 7. Limpieza
rm -f "$bam_raw" "$bam_raw.bai" "$metrics_file"
