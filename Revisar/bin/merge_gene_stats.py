#!/usr/bin/env python3
"""
merge_gene_stats.py
Fusiona varios <sample>.gene_stats.tsv en un TSV global ancho por GENE.
Cabecera por fichero (per-sample) esperada:
  GENE    %COVERED (sample)    # NON COVERED VARIANTS (sample)    NON COVERED VARIANTS (sample)

Salida:
  GENE    <SAMPLE>_%COVERED (sample)    <SAMPLE>_# NON COVERED VARIANTS (sample)    <SAMPLE>_NON COVERED VARIANTS (sample)    ...

Nota: no requiere pandas; usa csv/stdlib. Genes se ordenan alfabéticamente para estabilidad.
"""

import argparse, os, csv, sys
from collections import defaultdict, OrderedDict

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True, help="Lista de *.gene_stats.tsv por muestra")
    ap.add_argument("--output", required=True, help="Ruta de salida del TSV global")
    return ap.parse_args()

def sample_from_path(p):
    base = os.path.basename(p)
    if base.endswith(".gene_stats.tsv"):
        return base[:-len(".gene_stats.tsv")]
    # fallback: quita extensión .tsv si la hay
    if base.endswith(".tsv"):
        return base[:-4]
    return os.path.splitext(base)[0]

def read_gene_table(path):
    """
    Devuelve OrderedDict[GENE] = [pct_covered, n_nc, list_nc]
    """
    data = OrderedDict()
    with open(path, "rt", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr, None)
        if not header or len(header) < 4 or header[0] != "GENE":
            sys.exit(f"Cabecera inesperada en {path}. Se esperaba 4 columnas con 'GENE' en la primera.")
        for row in rdr:
            if not row or len(row) < 4:
                continue
            gene = row[0]
            pct  = row[1]
            n_nc = row[2]
            list_nc = row[3]
            data[gene] = [pct, n_nc, list_nc]
    return data

def main():
    a = parse_args()
    inputs = a.inputs
    if not inputs:
        sys.exit("No hay inputs.")

    # lee todas las tablas en memoria (paneles suelen ser pequeños)
    per_sample = []
    all_genes = set()
    for p in inputs:
        sample = sample_from_path(p)
        tbl = read_gene_table(p)
        per_sample.append((sample, tbl))
        all_genes.update(tbl.keys())

    # orden estable de genes
    genes_sorted = sorted(all_genes)

    # construye cabecera global
    header = ["GENE"]
    for sample, _ in per_sample:
        header += [
            f"{sample}_%COVERED (sample)",
            f"{sample}_# NON COVERED VARIANTS (sample)",
            f"{sample}_NON COVERED VARIANTS (sample)"
        ]

    # escribe salida
    with open(a.output, "wt", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(header)
        for gene in genes_sorted:
            row = [gene]
            for sample, tbl in per_sample:
                if gene in tbl:
                    pct, n_nc, list_nc = tbl[gene]
                else:
                    pct, n_nc, list_nc = ("0.00", "0", "-")
                row += [pct, n_nc, list_nc]
            w.writerow(row)

if __name__ == "__main__":
    main()
