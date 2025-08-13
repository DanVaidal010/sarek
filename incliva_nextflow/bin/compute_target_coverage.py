#!/usr/bin/env python3
"""
compute_target_coverage.py --mode gene
- Lee mosdepth *.per-base.bed.gz (chr, start, end, depth)
- Lee capture.bed (chr, start, end, [amplicon_id], [gene])
- Agrupa por GENE y calcula %COVERED (>= cutoff)
- Opcional: annotated.csv (con GENE_SYMBOL y POSITION 'chrN:pos') -> detecta variantes en bases < cutoff

Salida TSV:
GENE    %COVERED (sample)    # NON COVERED VARIANTS (sample)    NON COVERED VARIANTS (sample)
"""

import argparse, gzip, sys, csv
from collections import defaultdict, namedtuple

Interval = namedtuple("Interval", ["start", "end", "gene"])

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["gene"], required=True)
    ap.add_argument("--per-base", required=True, help="mosdepth *.per-base.bed.gz")
    ap.add_argument("--capture-bed", required=True, help="BED: chr start end [amplicon_id] [gene?]")
    ap.add_argument("--cutoff", type=int, default=250, help="Depth cutoff (>=) to consider covered")
    ap.add_argument("--annotated", help="annotated.csv (opcional). Debe contener columnas GENE_SYMBOL y POSITION (chr:pos)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--sample", required=True)
    return ap.parse_args()

def read_capture_bed(path):
    """
    Devuelve dict por cromosoma con lista de Interval(start,end,gene) ordenados por start.
    Columna gene:
      - si hay 5ª, se usa como gene
      - si no, si hay 4ª, se usa como gene
      - si no, 'NA'
    """
    by_chr = defaultdict(list)
    with open(path, "rt") as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"): continue
            toks = ln.rstrip("\n").split("\t")
            if len(toks) < 3: continue
            chrom = toks[0]
            try:
                start = int(toks[1]); end = int(toks[2])
            except ValueError:
                continue
            gene = "NA"
            if len(toks) >= 5 and toks[4]:
                gene = toks[4]
            elif len(toks) >= 4 and toks[3]:
                gene = toks[3]
            by_chr[chrom].append(Interval(start, end, gene))
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda iv: iv.start)
    return by_chr

def sweep_per_base(per_base_gz, capture_by_chr, cutoff, collect_lowcov_positions=False):
    """
    Recorre per-base y acumula por gen:
      totals[gene]: bases target totales
      covered[gene]: bases con depth >= cutoff
    Si collect_lowcov_positions=True, también devuelve:
      lowcov[gene][chrom] = set(positions<cutoff)
    """
    totals  = defaultdict(int)
    covered = defaultdict(int)
    lowcov  = defaultdict(lambda: defaultdict(set)) if collect_lowcov_positions else None

    with gzip.open(per_base_gz, "rt") as fh:
        current_chr = None
        intervals = []
        i = 0  # puntero al primer intervalo potencialmente solapante

        for ln in fh:
            if not ln.strip() or ln.startswith("#"): continue
            toks = ln.rstrip("\n").split("\t")
            if len(toks) < 4: continue
            chrom = toks[0]
            try:
                start = int(toks[1]); end = int(toks[2]); depth = float(toks[3])
            except ValueError:
                continue

            # mosdepth per-base usa [start,end) con longitud 1 base -> usamos 'pos=start'
            pos = start

            # cuando cambia de cromosoma, cargamos intervalos
            if chrom != current_chr:
                current_chr = chrom
                intervals = capture_by_chr.get(chrom, [])
                i = 0

            if not intervals:
                continue

            # avanza hasta el primer intervalo con end > pos
            while i < len(intervals) and intervals[i].end <= pos:
                i += 1

            # recorre intervalos solapantes (start <= pos < end)
            j = i
            while j < len(intervals) and intervals[j].start <= pos:
                iv = intervals[j]
                if pos < iv.end:
                    gene = iv.gene
                    totals[gene] += 1
                    if depth >= cutoff:
                        covered[gene] += 1
                    elif collect_lowcov_positions:
                        lowcov[gene][chrom].add(pos)
                j += 1

    return totals, covered, lowcov

def read_annotated_variants(path):
    """
    Lee annotated.csv y devuelve lista de (gene, chrom, pos, ident_str)
    Formato esperado (Excel de INCLIVA):
      - GENE_SYMBOL
      - POSITION con forma 'chrN:pos'
      - REF, ALT para formar un identificador opcional
    """
    variants = []
    with open(path, "rt", newline="") as fh:
        rdr = csv.DictReader(fh)
        for row in rdr:
            gene = row.get("GENE_SYMBOL") or row.get("Gene") or row.get("GENE") or "NA"
            pos_field = row.get("POSITION") or row.get("Pos") or row.get("POSITION (hg38)")
            if not pos_field or ":" not in pos_field:
                continue
            chrom, spos = pos_field.split(":", 1)
            try:
                pos = int(spos)
            except ValueError:
                continue
            ref = row.get("REF") or ""
            alt = row.get("ALT") or ""
            ident = f"{chrom}:{pos}"
            if ref and alt:
                ident += f":{ref}>{alt}"
            variants.append((gene, chrom, pos, ident))
    return variants

def write_gene_tsv(out_path, totals, covered, noncovered_by_gene):
    genes = sorted(set(totals.keys()) | set(covered.keys()) | set(noncovered_by_gene.keys()))
    with open(out_path, "wt") as out:
        out.write("GENE\t%COVERED (sample)\t# NON COVERED VARIANTS (sample)\tNON COVERED VARIANTS (sample)\n")
        for gene in genes:
            t = totals.get(gene, 0)
            c = covered.get(gene, 0)
            pct = 0.0 if t == 0 else 100.0 * c / t
            nc_list = noncovered_by_gene.get(gene, [])
            if nc_list:
                out.write(f"{gene}\t{pct:.2f}\t{len(nc_list)}\t{';'.join(nc_list)}\n")
            else:
                out.write(f"{gene}\t{pct:.2f}\t0\t-\n")

def main():
    a = parse_args()
    if a.mode != "gene":
        sys.exit("Only --mode gene is supported here.")

    capture_by_chr = read_capture_bed(a.capture_bed)

    # Si hay annotated.csv, recolectamos posiciones < cutoff para marcar variantes no cubiertas
    collect_lowcov = bool(a.annotated)
    totals, covered, lowcov = sweep_per_base(a.per_base, capture_by_chr, a.cutoff, collect_lowcov_positions=collect_lowcov)

    noncovered_by_gene = defaultdict(list)
    if a.annotated and lowcov is not None:
        variants = read_annotated_variants(a.annotated)
        # Por cada variante, si su pos pertenece a lowcov del gen/chrom, la marcamos como "non covered"
        for gene, chrom, pos, ident in variants:
            if pos in lowcov.get(gene, {}).get(chrom, set()):
                noncovered_by_gene[gene].append(ident)

    write_gene_tsv(a.output, totals, covered, noncovered_by_gene)

if __name__ == "__main__":
    main()
