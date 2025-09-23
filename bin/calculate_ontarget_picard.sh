#!/usr/bin/env bash
set -euo pipefail

# === INPUTS ===
cram_path="$1"
fastq_r1="$2"
output_dir="$3"
sample_id="$4"

# === CONFIG ===
reference_fasta="/Incliva/Sarek/Genomas.Ref/Homo_sapiens_assembly38.fasta"
PICARD_IMG="/Incliva/Sarek/Containers.Singularity/picard:3.4.0--hdfd78af_0"

# Binds robustos (libz suele estar en /lib/x86_64-linux-gnu)
BIND="/Incliva:/Incliva,/lib:/lib,/lib64:/lib64,/usr/lib:/usr/lib,/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu,/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu"

# Función wrapper para Picard con entorno saneado
run_picard() {
  SINGULARITYENV_LC_ALL=C.UTF-8 \
  SINGULARITYENV_LANG=C.UTF-8 \
  singularity exec -B "$BIND" "$PICARD_IMG" picard "$@"
}

echo "🧪 Procesando %OnTargetNoDupReads para muestra: $sample_id"

prefix="${output_dir}/${sample_id}"
bam_raw="${prefix}.mapped.bam"
metrics_file="${prefix}_ontarget_alignment_metrics.txt"
result_file="${prefix}.ontarget_result.txt"

# CRAM -> BAM
samtools view -b -T "$reference_fasta" "$cram_path" > "$bam_raw"
samtools index "$bam_raw"

# Picard
run_picard CollectAlignmentSummaryMetrics \
  I="$bam_raw" \
  O="$metrics_file" \
  R="$reference_fasta" \
  VALIDATION_STRINGENCY=LENIENT

ontarget_reads=$(grep -m1 "^PAIR" "$metrics_file" | awk '{print $6}' || echo "0")

# FASTQ estándar: 4 líneas por read
raw_reads=$(zcat "$fastq_r1" | wc -l | awk '{printf "%d", $1/2}')

if [[ -z "$ontarget_reads" || -z "$raw_reads" || "$raw_reads" -eq 0 ]]; then
  on_target_pct="0.00"
else
  on_target_pct=$(awk "BEGIN { printf \"%.2f\", ($ontarget_reads / $raw_reads) * 100 }")
fi

echo -e "%OnTargetNoDupReads\t$on_target_pct" > "$result_file"
echo "✅ Resultado guardado en $result_file"

rm -f "$bam_raw" "$bam_raw.bai" "$metrics_file"
