#!/usr/bin/env bash
set -euo pipefail

# Uso:
#   calculate_kit_specificity_picard.sh <cram> <outdir> <sample> <reference_fasta> <capture_bed>
CRAM="$1"
OUTDIR="$2"
SAMPLE="$3"
REF_FASTA="$4"
CAPTURE_BED="$5"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

export LC_ALL=C

log="${SAMPLE}.kit_specificity.debug.log"
echo "[INFO] sample=${SAMPLE}" > "$log"
echo "[INFO] cram=${CRAM}" >> "$log"
echo "[INFO] ref=${REF_FASTA}" >> "$log"
echo "[INFO] capture_bed=${CAPTURE_BED}" >> "$log"

# ---------- Archivos intermedios ----------
bam_all="${SAMPLE}.global.bam"
bam_target="${SAMPLE}.ontarget.bam"
metrics_all="${SAMPLE}.global.aln_metrics.txt"
metrics_target="${SAMPLE}.ontarget.aln_metrics.txt"
result_file="${SAMPLE}.kit_specificity_result.txt"

# ---------- CRAM -> BAM ----------
samtools view -b -T "$REF_FASTA" "$CRAM" > "$bam_all"
samtools index "$bam_all"

# ---------- Subconjunto on-target ----------
samtools view -b -L "$CAPTURE_BED" "$bam_all" > "$bam_target"
samtools index "$bam_target"

# ---------- Ejecutar Picard ----------
picard CollectAlignmentSummaryMetrics \
  I="${bam_all}" \
  O="${metrics_all}" \
  R="${REF_FASTA}" \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

picard CollectAlignmentSummaryMetrics \
  I="${bam_target}" \
  O="${metrics_target}" \
  R="${REF_FASTA}" \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

# ---------- Localizar índice PF_READS_ALIGNED ----------
get_pf_reads() {
  local file="$1"
  local idx
  idx=$(awk 'BEGIN{FS="\t"} $1=="CATEGORY"{for(i=1;i<=NF;i++){if($i=="PF_READS_ALIGNED"){print i; exit}}}' "$file")
  awk -v i="$idx" 'BEGIN{FS="\t"} $1=="PAIR"{print $(i); exit}' "$file"
}

mapped_reads=$(get_pf_reads "$metrics_all")
ontarget_reads=$(get_pf_reads "$metrics_target")
echo "[INFO] mapped_reads=${mapped_reads:-NA}" >> "$log"
echo "[INFO] ontarget_reads=${ontarget_reads:-NA}" >> "$log"

# ---------- % Kit specificity ----------
specificity="NA"
if [[ -n "${mapped_reads}" && -n "${ontarget_reads}" && "$mapped_reads" =~ ^[0-9]+$ && "$mapped_reads" -gt 0 ]]; then
  specificity=$(awk -v num="$ontarget_reads" -v den="$mapped_reads" 'BEGIN{printf "%.2f", 100*num/den}')
fi

# Clip 0..100
if [[ "$specificity" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
  specificity=$(awk -v x="$specificity" 'BEGIN{if(x<0)x=0; if(x>100)x=100; printf "%.2f", x}')
else
  specificity="NA"
fi
echo "[INFO] specificity=${specificity}" >> "$log"

# ---------- Output ----------
echo -e "Kit_specificity\t${specificity}" > "$result_file"
echo "✅ Resultado guardado en ${result_file}"
