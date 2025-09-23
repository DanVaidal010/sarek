#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import argparse
import traceback
import os

def parse_txt_to_dict(file_path):
    data = {}
    try:
        with open(file_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    key = parts[0].strip()
                    val = parts[1].strip()
                    data[key] = val
    except Exception as e:
        print(f"❌ Error al procesar {file_path}: {e}")
        traceback.print_exc()
    return data

# === ARGUMENTOS ===
parser = argparse.ArgumentParser(description="Fusiona estadísticas individuales en un único CSV")
parser.add_argument("--sample", required=True)
parser.add_argument("--quartile", required=True)
parser.add_argument("--mean", required=True)
parser.add_argument("--insert", required=True)
parser.add_argument("--numvars", required=True)
parser.add_argument("--basic", required=True)
parser.add_argument("--ontarget", required=True)
parser.add_argument("--kit", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

# === COMBINACIÓN DE ESTADÍSTICAS ===
stats = {
    "Sample": args.sample
}
stats.update(parse_txt_to_dict(args.basic))
stats.update(parse_txt_to_dict(args.quartile))
stats.update(parse_txt_to_dict(args.mean))
stats.update(parse_txt_to_dict(args.insert))

# Leer NumVars
try:
    with open(args.numvars) as f:
        stats["NumVars"] = f.readline().strip()
except Exception as e:
    print(f"❌ Error al leer NumVars: {e}")
    stats["NumVars"] = "NA"

# Leer %OnTargetNoDupReads
stats["%OnTargetNoDupReads"] = "NA"
if os.path.exists(args.ontarget):
    try:
        with open(args.ontarget) as f:
            for line in f:
                if line.startswith("%OnTargetNoDupReads"):
                    stats["%OnTargetNoDupReads"] = line.strip().split('\t')[1]
                    break
    except Exception as e:
        print(f"❌ Error al leer %OnTargetNoDupReads: {e}")
        traceback.print_exc()
else:
    print(f"⚠️ Archivo no encontrado: {args.ontarget}")

# Leer Kit_specificity
stats["Kit_specificity"] = "NA"
if os.path.exists(args.kit):
    try:
        with open(args.kit) as f:
            for line in f:
                if line.startswith("Kit_specificity"):
                    stats["Kit_specificity"] = line.strip().split('\t')[1]
                    break
    except Exception as e:
        print(f"❌ Error al leer Kit_specificity: {e}")
        traceback.print_exc()
else:
    print(f"⚠️ Archivo no encontrado: {args.kit}")

# === ORDEN DE COLUMNAS PARA EL CSV FINAL ===
columns = [
    "Sample", "Rawdata", "%LowQReads", "%MappedReads", "%Paired_reads", "%DuplicateReads",
    "%OnTargetNoDupReads", "Kit_specificity",
    "Cov_1stQ", "Cov_3rdQ", "Cov_Median",
    "Nt_100x", "Nt_250x", "Nt_500x",
    "Mean_insert", "Insert_SD", "NumVars"
]

# === ESCRIBIR ARCHIVO CSV ===
with open(args.output, 'w', newline='') as out:
    writer = csv.DictWriter(out, fieldnames=columns)
    writer.writeheader()
    writer.writerow(stats)

print(f"✅ CSV generado correctamente: {args.output}")
