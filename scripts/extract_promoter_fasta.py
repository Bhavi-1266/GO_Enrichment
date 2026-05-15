#!/usr/bin/env python3
"""Extract promoter DNA into FASTA with headers compatible with find_nrf1_genes.py.

Header format: >chr@promo_start-promo_end|GENE::chr:tss-tss_end(strand)

Coordinates are 1-based inclusive for pyfaidx Fasta.get_seq.
TSS is taken from BED start column (same convention as make_tss_bed.py output).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pyfaidx import Fasta


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Promoter FASTA from TSS BED + reference genome")
    p.add_argument("--bed", required=True, help="BED from make_tss_bed.py")
    p.add_argument("--genome", required=True, help="Reference FASTA (e.g. hg38.fa, may be .gz)")
    p.add_argument("--output", required=True, help="Output FASTA path")
    p.add_argument("--upstream", type=int, default=2000, help="Bases upstream of TSS")
    p.add_argument("--downstream", type=int, default=200, help="Bases downstream of TSS")
    return p.parse_args()


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def main() -> None:
    args = parse_args()
    bed_path = Path(args.bed)
    out_path = Path(args.output)
    fa = Fasta(str(args.genome), as_raw=True)

    n = 0
    skipped = 0

    with bed_path.open("r", encoding="utf-8") as fin, out_path.open("w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                skipped += 1
                continue
            chrom = cols[0]
            tss = int(cols[1])
            tss_end = int(cols[2])
            strand = cols[5]
            if strand == "+":
                p0 = tss - args.upstream
                p1 = tss + args.downstream
            else:
                p0 = tss - args.downstream
                p1 = tss + args.upstream
            if p0 < 1:
                p0 = 1
            try:
                chrom_len = len(fa[chrom])
            except KeyError:
                skipped += 1
                continue
            if p1 > chrom_len:
                p1 = chrom_len
            if p0 > p1:
                skipped += 1
                continue
            try:
                raw = str(fa.get_seq(chrom, p0, p1))
            except Exception:
                skipped += 1
                continue
            seq = raw.upper()
            if strand == "-":
                seq = revcomp(seq)
            name_field = cols[3]
            parts = name_field.split("|", 1)
            if len(parts) < 2:
                gene = parts[0]
            else:
                gene = parts[1]
            hdr = f">{chrom}@{p0}-{p1}|{gene}::{chrom}:{tss}-{tss_end}({strand})"
            fout.write(hdr + "\n")
            for i in range(0, len(seq), 60):
                fout.write(seq[i : i + 60] + "\n")
            n += 1

    print(f"Wrote {n} promoter records to {out_path} (skipped {skipped}).")


if __name__ == "__main__":
    main()
