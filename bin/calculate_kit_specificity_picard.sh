#!/usr/bin/env bash
set -euo pipefail

# USO:
#   calculate_kit_specificity_picard.sh <cram> <outdir> <sample> <reference_fasta> <bed>
CRAM="$1"
OUTDIR="$2"
SAMPLE="$3"
REF_FASTA="$4"
BED="$5"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

if [[ ! -s "${BED}" ]]; then
  echo -e "Kit_specificity\tNA" > "${SAMPLE}.kit_specificity_result.txt"
  exit 0
fi

# Construir BAM global y on-target (CRAM necesita -T)
samtools view -T "${REF_FASTA}" -b "${CRAM}" > "${SAMPLE}.global.bam"
samtools view -T "${REF_FASTA}" -b -L "${BED}" "${CRAM}" > "${SAMPLE}.ontarget.bam"

# Picard para ambos
picard CollectAlignmentSummaryMetrics \
  I="${SAMPLE}.global.bam" \
  O="${SAMPLE}.global.aln_metrics.txt" \
  REFERENCE_SEQUENCE="${REF_FASTA}" \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

picard CollectAlignmentSummaryMetrics \
  I="${SAMPLE}.ontarget.bam" \
  O="${SAMPLE}.ontarget.aln_metrics.txt" \
  REFERENCE_SEQUENCE="${REF_FASTA}" \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

# PF_READS_ALIGNED (PAIR o FIRST+SECOND)
get_pf() {
  local f="$1"
  local x
  x=$(awk 'BEGIN{FS="\t"} $1=="PAIR"{print $10}' "$f" 2>/dev/null || true)
  if [[ -z "$x" || "$x" == "0" ]]; then
    local a b
    a=$(awk 'BEGIN{FS="\t"} $1=="FIRST_OF_PAIR"{print $10}'  "$f" 2>/dev/null || echo 0)
    b=$(awk 'BEGIN{FS="\t"} $1=="SECOND_OF_PAIR"{print $10}' "$f" 2>/dev/null || echo 0)
    x=$(awk -v aa="$a" -v bb="$b" 'BEGIN{print (aa+bb)}')
  fi
  echo "$x"
}

pf_global=$(get_pf "${SAMPLE}.global.aln_metrics.txt")
pf_target=$(get_pf "${SAMPLE}.ontarget.aln_metrics.txt")

pct="NA"
if [[ "$pf_global" -gt 0 ]]; then
  pct=$(awk -v a="$pf_target" -v g="$pf_global" 'BEGIN{printf "%.2f", 100*a/g}')
fi

# Clip 0..100
if [[ "${pct}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
  pct=$(awk -v x="$pct" 'BEGIN{if(x<0)x=0; if(x>100)x=100; printf "%.2f", x}')
else
  pct="NA"
fi

echo -e "Kit_specificity\t${pct}" > "${SAMPLE}.kit_specificity_result.txt"
