#!/usr/bin/env bash
set -euo pipefail

# Uso:
# ./count_numvars.sh input.vcf output.txt
# Si output_file es "-", se imprime el número en stdout.

input_file="$1"     # Archivo VCF de entrada (con cabecera y variantes)
output_file="$2"    # Archivo de salida o "-" para imprimir por pantalla

# Contar el número de variantes:
#  - grep -v "^#" : excluye las líneas de cabecera del VCF
#  - wc -l        : cuenta el número de líneas (variantes)
num_vars=$(grep -v "^#" "$input_file" | wc -l)

# Si la salida es "-", mostrar por stdout
if [[ "$output_file" == "-" ]]; then
  echo "$num_vars"
else
  # Si no, guardar en archivo
  echo "$num_vars" > "$output_file"
fi
