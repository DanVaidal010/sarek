#!/usr/bin/env bash
set -euo pipefail
# USO:
#   count_numvars.sh <bcftools_stats_or_empty> [vcf_gz_optional]
STATS="$1"
VCF="$2"

if [[ -s "$STATS" ]]; then
  awk '/^SN/ && /number of records:/ {print $NF; exit}' "$STATS"
  exit 0
fi

if [[ -s "$VCF" ]]; then
  bcftools view -H "$VCF" 2>/dev/null | wc -l | awk '{print $1}'
  exit 0
fi

echo "NA"
