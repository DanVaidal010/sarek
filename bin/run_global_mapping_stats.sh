#!/bin/bash
set -euo pipefail

# === COMPROBAR ARGUMENTO ===
if [[ $# -lt 1 ]]; then
  echo "❌ Debes proporcionar el nombre de la muestra como argumento."
  echo "➡️ Uso: ./run_global_mapping_stats.sh test_sample_1"
  exit 1
fi

sample="$1"

# === PARÁMETROS ===
script_dir="/Incliva/Sarek/pid/global_mapping_stats/bin"
output_base="/Incliva/Sarek/pid/global_mapping_stats/output"
output_dir="${output_base}/${sample}"
json_file="/Incliva/Sarek/Resultados/reports/fastp/${sample}/${sample}-lane1.fastp.json"
cram_stats="/Incliva/Sarek/Resultados/reports/samtools/${sample}/${sample}.recal.cram.stats"
thresholds_bed="/Incliva/Sarek/Resultados/reports/mosdepth/${sample}/${sample}.sorted.thresholds.bed.gz"
perbase_bed="/Incliva/Sarek/Resultados/reports/mosdepth/${sample}/${sample}.recal.per-base.bed.gz"
bcftools_stats="/Incliva/Sarek/Resultados/reports/bcftools/mutect2/${sample}/${sample}.mutect2.filtered.bcftools_stats.txt"
panel_bed="/Incliva/targetseq-master/targetseq-master/panel_designs/TOX4/capture.bed"
sample_cram="/Incliva/Sarek/Resultados/preprocessing/recalibrated/${sample}/${sample}.recal.cram"
reference="/Incliva/Sarek/Genomas.Ref/Homo_sapiens_assembly38.fasta"
picard_img="/Incliva/Sarek/Containers.Singularity/picard:3.4.0--hdfd78af_0"

# === MANUAL OVERRIDE PARA TESTING ===
if [[ "$sample" == "test_sample_1" ]]; then
  fastq_r1="/Incliva/Sarek/Muestras/TOX4/TOX4-DDPD-1648_S7_L001_R1_001.fastq.gz"
else
  fastq_r1=$(find /Incliva/Sarek/Muestras/TOX4/ -type f -name "*R1*.fastq.gz" | grep -i "$sample" | head -n 1)
fi

# === Verificar existencia del archivo FASTQ R1 ===
if [[ ! -f "$fastq_r1" ]]; then
  echo "❌ FASTQ R1 no encontrado: $fastq_r1"
  exit 1
else
  echo "✅ FASTQ R1 detectado: $fastq_r1"
fi

# === Crear directorio de salida ===
mkdir -p "$output_dir"

# === Definir nombres de archivo de salida ===
basic_output="${output_dir}/${sample}.fastq_summary.txt"
cuartil_output="${output_dir}/${sample}.coverage_quartiles.txt"
mean_output="${output_dir}/${sample}.coverage_thresholds.txt"
insert_output="${output_dir}/${sample}.insert_metrics.txt"
numvars_output="${output_dir}/${sample}.variant_count.txt"
ontarget_output="${output_dir}/${sample}.ontarget_result.txt"
kit_output="${output_dir}/${sample}.kit_specificity_result.txt"
final_csv="${output_dir}/${sample}_mapping_stats.csv"

# === INICIO DEL PIPELINE ===
echo "▶️ [1/8] Estadísticas básicas (fastp + samtools)..."
bash "$script_dir/summarize_fastq_cram_stats.sh" "$json_file" "$cram_stats" "$sample" "$basic_output"

echo "▶️ [2/8] Cuartiles de cobertura..."
python3 "$script_dir/extract_coverage_quartiles.py" -p "$perbase_bed" -c "$panel_bed" -o "$cuartil_output"

echo "▶️ [3/8] Coberturas 100x/250x/500x..."
python3 "$script_dir/coverage_depth_thresholds_100_250_500.py" -i "$thresholds_bed" -o "$mean_output"

echo "▶️ [4/8] Inserto medio y desviación estándar..."
python3 "$script_dir/calculate_insert_metrics.py" -j "$json_file" -o "$insert_output"

echo "▶️ [5/8] Conteo de variantes..."
bash "$script_dir/count_filtered_variants.sh" "$bcftools_stats" "$numvars_output"

echo "▶️ [6/8] %OnTargetNoDupReads con Picard..."
bash "$script_dir/calculate_ontarget_picard.sh" \
  "$sample_cram" \
  "$fastq_r1" \
  "$output_dir" \
  "$sample" \
  "$reference" \
  "$picard_img"

echo "▶️ [7/8] Kit Specificity con Picard..."
bash "$script_dir/calculate_kit_specificity_picard.sh" \
  "$sample_cram" \
  "$output_dir" \
  "$sample" \
  "$reference" \
  "$panel_bed" \
  "$picard_img"


echo "📦 [8/8] Generando CSV final..."
python3 "$script_dir/combine_mapping_statistics.py" \
  --sample "$sample" \
  --quartile "$cuartil_output" \
  --mean "$mean_output" \
  --insert "$insert_output" \
  --numvars "$numvars_output" \
  --basic "$basic_output" \
  --ontarget "$ontarget_output" \
  --kit "$kit_output" \
  --output "$final_csv"

echo "✅ CSV generado en: $final_csv"
