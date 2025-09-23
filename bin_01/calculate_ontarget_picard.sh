#!/usr/bin/env bash
set -euo pipefail

# Uso:
#   calculate_ontarget_picard.sh <cram> <fastq_r1_ignored> <outdir> <sample> <reference_fasta> <cram_stats>
CRAM="$1"
R1_IGNORED="$2"
OUTDIR="$3"
SAMPLE="$4"
REF_FASTA="$5"
CRAM_STATS="$6"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

export LC_ALL=C

log="${SAMPLE}.ontarget.debug.log"
echo "[INFO] sample=${SAMPLE}" > "$log"
echo "[INFO] cram=${CRAM}" >> "$log"
echo "[INFO] ref=${REF_FASTA}" >> "$log"
echo "[INFO] cram_stats=${CRAM_STATS}" >> "$log"

# ---------- Denominador: raw total sequences ----------
raw_total=$(awk -F'\t' '/^SN\traw total sequences:/ {gsub(/^[ \t]+/,"",$3); print $3; exit}' "${CRAM_STATS}" 2>/dev/null)
if [[ -z "$raw_total" ]]; then
  raw_total=$(awk -F':' '/raw total sequences:/ {gsub(/^[ \t]+/,"",$2); print $2; exit}' "${CRAM_STATS}" 2>/dev/null || echo "NA")
fi
echo "[INFO] raw_total=${raw_total}" >> "$log"

# ---------- Ejecutar Picard ----------
picard CollectAlignmentSummaryMetrics \
  I="${CRAM}" \
  O="${SAMPLE}.picard.aln_metrics.txt" \
  REFERENCE_SEQUENCE="${REF_FASTA}" \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

if [[ ! -s "${SAMPLE}.picard.aln_metrics.txt" ]]; then
  echo "[ERROR] Picard no generó ${SAMPLE}.picard.aln_metrics.txt" >> "$log"
  echo -e "%OnTargetNoDupReads\tNA" > "${SAMPLE}.ontarget_result.txt"
  exit 0
fi

# ---------- Localizar índice PF_READS_ALIGNED ----------
col_pf_aligned=$(awk 'BEGIN{FS="\t"} $1=="CATEGORY"{for(i=1;i<=NF;i++){if($i=="PF_READS_ALIGNED"){print i; exit}}}' "${SAMPLE}.picard.aln_metrics.txt")
echo "[INFO] col_pf_aligned=${col_pf_aligned:-NA}" >> "$log"

# ---------- Tomar PAIR; si no, FIRST+SECOND ----------
get_val_by_cat() {
  local cat="$1" idx="$2"
  awk -v c="$cat" -v i="$idx" 'BEGIN{FS="\t"} $1==c{print $(i); exit}' "${SAMPLE}.picard.aln_metrics.txt"
}

pf_aligned="$(get_val_by_cat PAIR "${col_pf_aligned}")"
if [[ -z "$pf_aligned" || "$pf_aligned" == "0" ]]; then
  first="$(get_val_by_cat FIRST_OF_PAIR "${col_pf_aligned}")"
  second="$(get_val_by_cat SECOND_OF_PAIR "${col_pf_aligned}")"
  first="${first:-0}"; second="${second:-0}"
  pf_aligned=$(awk -v a="$first" -v b="$second" 'BEGIN{print a+b}')
fi
echo "[INFO] pf_aligned=${pf_aligned:-NA}" >> "$log"

# ---------- %OnTargetNoDupReads ----------
pct="NA"
if [[ "$raw_total" != "NA" && "$raw_total" =~ ^[0-9]+$ && "$raw_total" -gt 0 ]]; then
  pct=$(awk -v num="$pf_aligned" -v den="$raw_total" 'BEGIN{printf "%.2f", 100*num/den}')
fi

# Clip 0..100
if [[ "$pct" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
  pct=$(awk -v x="$pct" 'BEGIN{if(x<0)x=0; if(x>100)x=100; printf "%.2f", x}')
else
  pct="NA"
fi
echo "[INFO] pct=${pct}" >> "$log"

echo -e "%OnTargetNoDupReads\t${pct}" > "${SAMPLE}.ontarget_result.txt"
echo "✅ Resultado guardado en ${SAMPLE}.ontarget_result.txt"
