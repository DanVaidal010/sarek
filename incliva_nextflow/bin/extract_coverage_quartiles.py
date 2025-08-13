#!/usr/bin/env python3 
import pybedtools
import pandas as pd
import numpy as np
import argparse

def obtener_cobertura_intersect(perbase_path, capture_bed_path):
    """
    Intersecta el archivo de cobertura por base (`perbase.bed.gz`)
    con las regiones del panel de captura (`capture.bed`),
    expandiendo los intervalos a posiciones individuales y devolviendo
    la cobertura por posición.
    """
    # Cargar archivos como objetos BedTool
    perbase = pybedtools.BedTool(perbase_path)
    capture = pybedtools.BedTool(capture_bed_path)

    # Intersectar para quedarnos solo con las posiciones dentro del panel
    intersect = perbase.intersect(capture, wa=True)

    # Pasar el resultado a DataFrame con nombres de columnas
    intersect_df = pd.read_csv(intersect.fn, sep='\t', header=None, names=['chrom','start','end','coverage'])

    # Eliminar duplicados exactos
    intersect_df = intersect_df.drop_duplicates()

    # Desdoblar intervalos en posiciones individuales
    intersect_df['pos'] = intersect_df.apply(lambda row: list(range(row['start'], row['end'])), axis=1)
    expanded_df = intersect_df.explode('pos').reset_index(drop=True)

    # Guardar archivo expandido (opcional, útil para depuración)
    expanded_df[['chrom', 'pos', 'pos', 'coverage']].to_csv('intersect_expanded.bed', sep='\t', header=False, index=False)

    # Devolver la serie de coberturas por posición
    return expanded_df['coverage']

def calcular_cuartiles(perbase_path, capture_bed_path):
    """
    Calcula el primer cuartil (Q1), tercer cuartil (Q3) y la mediana
    de la cobertura en las regiones del panel.
    """
    coverage = obtener_cobertura_intersect(perbase_path, capture_bed_path)
    if coverage.empty:
        raise ValueError("No se encontraron posiciones válidas después de la intersección.")
    q1 = np.percentile(coverage, 25)
    q3 = np.percentile(coverage, 75)
    median = np.median(coverage)
    return q1, q3, median

if __name__ == "__main__":
    # Parser de argumentos CLI
    parser = argparse.ArgumentParser(description="Cuartiles de cobertura con pybedtools")
    parser.add_argument("-p", "--perbase", required=True, help="Archivo per-base.bed.gz")
    parser.add_argument("-c", "--capture", required=True, help="Archivo capture.bed")
    parser.add_argument("-o", "--output", required=True, help="Archivo salida")
    args = parser.parse_args()

    # Calcular cuartiles
    q1, q3, med = calcular_cuartiles(args.perbase, args.capture)

    # Guardar resultados en archivo de texto
    with open(args.output, 'w') as f:
        f.write(f"Cov_1stQ\t{round(q1,2)}\n")
        f.write(f"Cov_3rdQ\t{round(q3,2)}\n")
        f.write(f"Cov_Median\t{round(med,2)}\n")

    print(f"✅ Cuartiles guardados en {args.output}")
