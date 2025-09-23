#!/bin/bash
set -euo pipefail

# === INPUTS ===
cram_path="$1"
fastq_r1="$2"
output_dir="$3"
sample_id="$4"

# === CONFIGURACIÓN ===
reference_fasta="/Incliva/Sarek/Genomas.Ref/Homo_sapiens_assembly38.fasta"

echo "🧪 Procesando %OnTargetNoDupReads para muestra: $sample_id (modo DEDUPLICADO)"

# === RUTAS DE SALIDA ===
prefix="${output_dir}/${sample_id}"
bam_raw="${prefix}.mapped.bam"
metrics_file="${prefix}_ontarget_alignment_metrics.txt"
result_file="${prefix}.ontarget_result.txt"

# === Convertir CRAM a BAM ===
samtools view -b -T "$reference_fasta" "$cram_path" > "$bam_raw"
samtools index "$bam_raw"

# === Extraer métricas de BAM ===
picard CollectAlignmentSummaryMetrics \
    I="$bam_raw" \
    O="$metrics_file" \
    R="$reference_fasta" \
    VALIDATION_STRINGENCY=LENIENT

# === Obtener número de lecturas no duplicadas (on-target) ===
ontarget_reads=$(grep -m 1 "^PAIR" "$metrics_file" | awk '{print $6}' || echo "0")

# === Calcular número total de raw reads desde FASTQ R1 ===
raw_reads=$(zcat "$fastq_r1" | echo $((`wc -l` / 2)))

# === Calcular porcentaje ===
if [[ -z "$ontarget_reads" || -z "$raw_reads" || "$raw_reads" -eq 0 ]]; then
  on_target_pct="0.00"
else
  on_target_pct=$(awk "BEGIN { printf \"%.2f\", ($ontarget_reads / $raw_reads) * 100 }")
fi

# === Guardar resultado ===
echo -e "%OnTargetNoDupReads\t$on_target_pct" > "$result_file"
echo "✅ Resultado guardado en $result_file"

# === Limpieza opcional ===
rm -f "$bam_raw" "$bam_raw.bai" "$metrics_file"
