#!/usr/bin/env bash
set -euo pipefail
# USO: count_filtered_variants.sh <bcftools_stats_or_empty> [vcf_gz_optional]

STATS="$1"
VCF="${2:-}"

num="NA"

if [[ -s "$STATS" ]]; then
  num=$(awk '/^SN/ && /number of records:/ {print $NF; exit}' "$STATS")
elif [[ -n "$VCF" && -s "$VCF" ]]; then
  num=$(bcftools view -H "$VCF" 2>/dev/null | wc -l)
fi

echo -e "NumVars\t${num}"
