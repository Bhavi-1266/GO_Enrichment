#!/usr/bin/env python3
"""Extract unique gene symbols whose promoter sequence contains NRF1 motif.

Motif from class note: GCGC..GCGC, where '.' can be any nucleotide.
Regex used: GCGC[ACGT]{2}GCGC
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


MOTIF = re.compile(r"GCGC[ACGT]{2}GCGC", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Find genes containing NRF1 motif in FASTA")
    parser.add_argument("--input", required=True, help="Input FASTA of promoter sequences")
    parser.add_argument("--output", required=True, help="Output gene list (one symbol per line)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    in_path = Path(args.input)
    out_path = Path(args.output)

    genes: set[str] = set()
    current_gene: str | None = None
    total_hits = 0

    with in_path.open("r", encoding="utf-8") as fin:
        for raw in fin:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Header expected format: >chr@start-end|GENE::chr:start-end(+/-)
                try:
                    current_gene = line.split("|", 1)[1].split("::", 1)[0]
                except IndexError:
                    current_gene = None
            else:
                if current_gene and MOTIF.search(line):
                    genes.add(current_gene)
                    total_hits += 1

    with out_path.open("w", encoding="utf-8") as fout:
        for gene in sorted(genes):
            fout.write(gene + "\n")

    print(f"Unique genes with motif: {len(genes)}")
    print(f"Total motif hits: {total_hits}")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
