#!/usr/bin/env python3
"""Convert human gene annotation TSV.GZ to TSS BED format for promoter analysis.

Expected input columns (by position, from assignment docs/examples):
- col5: chromosome (e.g., 1, X, MT)
- col6: strand code (1 or -1)
- col7: gene symbol
- col8: transcription start site (TSS)
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create TSS BED from human annotation TSV.GZ")
    parser.add_argument("--input", required=True, help="Path to human_gene_annotation.tsv.gz")
    parser.add_argument("--output", required=True, help="Path to output BED (genes_tss_clean.bed)")
    return parser.parse_args()


def to_chrom(raw: str) -> str:
    chrom = f"chr{raw}"
    if chrom == "chrMT":
        return "chrM"
    return chrom


def main() -> None:
    args = parse_args()
    in_path = Path(args.input)
    out_path = Path(args.output)

    written = 0
    skipped = 0

    with gzip.open(in_path, "rt") as fin, out_path.open("w", encoding="utf-8") as fout:
        _ = fin.readline()  # header
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                skipped += 1
                continue

            chrom_raw = cols[4].strip()
            strand_code = cols[5].strip()
            gene = cols[6].strip()
            tss_raw = cols[7].strip()

            if not chrom_raw or not gene or not tss_raw:
                skipped += 1
                continue

            try:
                tss = int(float(tss_raw))
            except ValueError:
                skipped += 1
                continue

            strand = "+" if strand_code in {"1", "+"} else "-"
            chrom = to_chrom(chrom_raw)
            tss_end = tss + 1
            name = f"{chrom}@{tss}-{tss_end}|{gene}"
            fout.write(f"{chrom}\t{tss}\t{tss_end}\t{name}\t.\t{strand}\n")
            written += 1

    print(f"Done. Written: {written}, Skipped: {skipped}")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
