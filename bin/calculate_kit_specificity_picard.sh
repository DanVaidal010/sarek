#!/usr/bin/env bash
set -euo pipefail

# === INPUTS ===
cram_path="$1"        # CRAM recalibrado
reference_fasta="$2"  # FASTA de referencia
capture_bed="$3"      # BED del panel de captura
sample_id="$4"        # ID de muestra
output_dir="$5"       # Carpeta de salida

# === OUTPUTS ===
prefix="${output_dir}/${sample_id}"
bam="${prefix}.bam"
bam_target="${prefix}.ontarget.bam"
metrics_all="${prefix}_mapped_alignment_metrics.txt"
metrics_target="${prefix}_ontarget_alignment_metrics.txt"
result_file="${prefix}.kit_specificity_result.txt"

echo "🧪 Procesando Kit_specificity para muestra: $sample_id"

# 1. CRAM → BAM global
samtools view -b -T "$reference_fasta" "$cram_path" > "$bam"
samtools index "$bam"

# 2. BAM filtrado por panel (on-target)
samtools view -b -L "$capture_bed" "$bam" > "$bam_target"
samtools index "$bam_target"

# 3. Picard sobre BAM global
picard CollectAlignmentSummaryMetrics \
    I="$bam" \
    O="$metrics_all" \
    R="$reference_fasta" \
    VALIDATION_STRINGENCY=LENIENT

# 4. Picard sobre BAM on-target
picard CollectAlignmentSummaryMetrics \
    I="$bam_target" \
    O="$metrics_target" \
    R="$reference_fasta" \
    VALIDATION_STRINGENCY=LENIENT

# 5. Extraer lecturas
mapped_reads=$(grep -m1 "^PAIR" "$metrics_all" | awk '{print $6}' || echo "0")
ontarget_reads=$(grep -m1 "^PAIR" "$metrics_target" | awk '{print $6}' || echo "0")

# 6. Calcular especificidad
if [[ -z "$mapped_reads" || -z "$ontarget_reads" || "$mapped_reads" -eq 0 ]]; then
    specificity="0.00"
else
    specificity=$(awk "BEGIN { printf \"%.2f\", ($ontarget_reads / $mapped_reads) * 100 }")
fi

# 7. Guardar resultado
echo -e "Kit_specificity\t$specificity" > "$result_file"
echo "✅ Resultado guardado en $result_file"

# 8. Limpieza
rm -f "$bam" "$bam.bai" "$bam_target" "$bam_target.bai" "$metrics_all" "$metrics_target"
