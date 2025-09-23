#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import traceback
import gzip

def calcular_porcentajes(path):
    """
    Calcula el % de bases cubiertas a >=100x, >=250x y >=500x
    a partir de un archivo .thresholds.bed.gz.

    El archivo debe contener:
    chr, start, end, cov_100, cov_250, cov_500
    """
    try:
        total_bases = 0
        bases_100x = 0
        bases_250x = 0
        bases_500x = 0

        with gzip.open(path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue  # Saltar cabecera
                parts = line.strip().split('\t')

                # Calcular longitud del segmento
                length = int(parts[2]) - int(parts[1])

                # Columnas con nº de bases que cumplen umbral
                cov_100 = int(parts[4])
                cov_250 = int(parts[5])
                cov_500 = int(parts[6])

                total_bases += length
                bases_100x += cov_100
                bases_250x += cov_250
                bases_500x += cov_500

        if total_bases == 0:
            raise ValueError("No se pudo calcular cobertura: total_bases = 0")

        # Calcular porcentajes
        return {
            "Nt_100x": round((bases_100x / total_bases) * 100, 2),
            "Nt_250x": round((bases_250x / total_bases) * 100, 2),
            "Nt_500x": round((bases_500x / total_bases) * 100, 2),
        }

    except Exception as e:
        print(f"❌ Error al calcular coberturas: {e}")
        traceback.print_exc()
        return {
            "Nt_100x": None,
            "Nt_250x": None,
            "Nt_500x": None,
        }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calcula % de cobertura a 100x, 250x y 500x.")
    parser.add_argument("-i", "--input", required=True, help="Ruta al archivo .thresholds.bed.gz")
    parser.add_argument("-s", "--sample", required=False, help="Nombre de la muestra (opcional)")
    parser.add_argument("-o", "--output", required=True, help="Archivo de salida")
    args = parser.parse_args()

    resultados = calcular_porcentajes(args.input)

    with open(args.output, 'w') as out:
        if args.sample:
            out.write(f"Sample\t{args.sample}\n")
        for clave, valor in resultados.items():
            out.write(f"{clave}\t{valor}\n")

    print(f"✅ Coberturas guardadas correctamente en {args.output}")
