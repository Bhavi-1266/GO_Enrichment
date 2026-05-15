#!/usr/bin/env python3
"""Build human_gene_annotation.tsv.gz from a GENCODE GTF (gene rows).

Output columns (tab-separated, gzip):
  col1..col4: placeholders (match assignment layout)
  chrom: chromosome without 'chr' prefix (1, X, MT) for make_tss_bed.py
  strand: 1 or -1
  gene_symbol: gene_name from GTF
  tss: 1-based transcription start (start for + strand, end for - strand)
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="GENCODE GTF → human_gene_annotation.tsv.gz")
    p.add_argument("--input", required=True, help="Path to .gtf or .gtf.gz")
    p.add_argument("--output", required=True, help="Path to output .tsv.gz")
    p.add_argument(
        "--types",
        default="protein_coding,lncRNA",
        help="Comma-separated gene_type values to keep (empty = all genes with gene_name)",
    )
    return p.parse_args()


def chrom_field(seqname: str) -> str:
    s = seqname.strip()
    if s.startswith("chr"):
        s = s[3:]
    if s == "M":
        return "M"
    return s


def parse_gene_types(raw: str) -> set[str] | None:
    raw = raw.strip()
    if not raw:
        return None
    return {t.strip() for t in raw.split(",") if t.strip()}


def main() -> None:
    args = parse_args()
    in_path = Path(args.input)
    out_path = Path(args.output)
    allow_types = parse_gene_types(args.types)

    def open_gtf() -> gzip.GzipFile | object:
        if str(in_path).endswith(".gz"):
            return gzip.open(in_path, "rt", encoding="utf-8", errors="replace")
        return in_path.open("r", encoding="utf-8", errors="replace")

    written = 0
    skipped = 0

    with open_gtf() as fin, gzip.open(
        out_path, "wt", encoding="utf-8", newline="\n"
    ) as fout:
        fout.write("c1\tc2\tc3\tc4\tchrom\tstrand\tgene_symbol\ttss\n")
        for line in fin:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            chrom = chrom_field(cols[0])
            start, end = int(cols[3]), int(cols[4])
            strand_sym = cols[6]
            attrs = cols[8]
            m = GENE_NAME_RE.search(attrs)
            if not m:
                skipped += 1
                continue
            gene = m.group(1)
            gtype_m = re.search(r'gene_type "([^"]+)"', attrs)
            gtype = gtype_m.group(1) if gtype_m else ""
            if allow_types is not None and gtype not in allow_types:
                skipped += 1
                continue
            strand_code = "1" if strand_sym == "+" else "-1"
            tss = start if strand_sym == "+" else end
            fout.write(f".\t.\t.\t.\t{chrom}\t{strand_code}\t{gene}\t{tss}\n")
            written += 1

    print(f"Wrote {written} gene rows to {out_path} (skipped {skipped}).")


if __name__ == "__main__":
    main()
