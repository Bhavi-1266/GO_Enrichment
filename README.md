# GO enrichment from NRF1 promoter motif scan

This repository runs a small end-to-end workflow: build human transcription start sites (TSS) from GENCODE, extract promoter sequences from hg38, find genes whose promoters match an NRF1-like motif (`GCGCNNGCGC`), then run Gene Ontology biological process (GO BP) enrichment on those genes using R/Bioconductor.

## What it does

1. **Annotation** — Downloads GENCODE basic GTF (default: release 45) and builds `input/human_gene_annotation.tsv.gz` (protein-coding genes).
2. **TSS BED** — Writes `output/genes_tss_clean.bed` for promoter extraction.
3. **Reference** — Downloads UCSC hg38 FASTA (large, cached under `input/`) and decompresses it to a plain `.fa` file for `pyfaidx`.
4. **Promoters** — Extracts a ±2 kb / +200 bp window around each TSS into `output/promoters.fa` (uses [pyfaidx](https://github.com/mdshw5/pyfaidx)).
5. **Motif scan** — Lists unique gene symbols with the motif in `output/nrf1_genes.txt`.
6. **GO enrichment** — Maps symbols to Entrez IDs and runs `clusterProfiler::enrichGO` (BP, BH-adjusted *p* & *q* ≤ 0.05). Writes `output/go_results.csv` and PNG plots when terms are found.

Orchestration lives in `scripts/run_pipeline.sh`.

## Requirements

- **Bash**, **wget** (or adjust downloads manually)
- **Python 3** with dependencies from `requirements.txt` (see below)
- **BioPython** (installed via `requirements.txt`) is required because `pyfaidx` reads the gzipped hg38 FASTA
- **Docker** — GO enrichment runs in `ghcr.io/bioconductor/bioconductor_docker:RELEASE_3_19` so you do not need a local R/Bioconductor install for that step

Optional: **R** with `clusterProfiler` and `org.Hs.eg.db` if you prefer to run `go_enrichment.R` on the host instead of in Docker.

## Quick start

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
bash scripts/run_pipeline.sh
```

If you already created `.venv`, just run:

```bash
source .venv/bin/activate
pip install -r requirements.txt
bash scripts/run_pipeline.sh
```

The first run downloads GENCODE GTF and hg38 (several GB for the genome) and then decompresses `input/hg38.fa.gz` into `input/hg38.fa`. Subsequent runs reuse files under `input/` and `output/`.

## Outputs (`output/`)

| File | Description |
|------|-------------|
| `genes_tss_clean.bed` | TSS intervals for promoter extraction |
| `promoters.fa` | Promoter sequences per gene |
| `nrf1_genes.txt` | One HGNC symbol per line (motif-positive promoters) |
| `go_results.csv` | GO BP enrichment table from clusterProfiler |
| `go_dotplot.png`, `go_barplot.png`, `go_cnetplot.png` | Plots (only if enrichment returns terms) |

## Configuration

Environment variables read by `scripts/run_pipeline.sh`:

| Variable | Default | Purpose |
|----------|---------|---------|
| `GTF_URL` | GENCODE v45 basic GTF.gz | Override GTF source |
| `HG38_URL` | UCSC hg38 `hg38.fa.gz` | Override reference FASTA |

The Python stage uses `.venv/bin/python3` when present; otherwise it falls back to `python3` on your `PATH`.

## Running GO enrichment alone

`go_enrichment.R` accepts optional CLI arguments: input gene list and output directory. Defaults match this pipeline.

```bash
Rscript go_enrichment.R path/to/genes.txt path/to/outdir
```

Inside Docker (same image as the shell script):

```bash
docker run --rm -v "$PWD:/work" -w /work ghcr.io/bioconductor/bioconductor_docker:RELEASE_3_19 \
  bash -lc 'Rscript -e "BiocManager::install(c(\"clusterProfiler\",\"org.Hs.eg.db\"), ask=FALSE, update=FALSE)" && Rscript go_enrichment.R output/nrf1_genes.txt output'
```

Human symbols are mapped with `org.Hs.eg.db`; enrichment is **GO BP** only, as implemented in `go_enrichment.R`.

## Scripts (reference)

| Script | Role |
|--------|------|
| `scripts/gtf_to_human_gene_annotation.py` | GTF → gzipped TSV for TSS extraction |
| `scripts/make_tss_bed.py` | Annotation TSV → BED |
| `scripts/extract_promoter_fasta.py` | BED + genome FASTA → promoter FASTA |
| `scripts/find_nrf1_genes.py` | Motif scan → gene symbol list |
| `scripts/run_pipeline.sh` | Full pipeline including Docker GO step |

## Notes

- Motif definition follows the in-script comment: `GCGC[ACGT]{2}GCGC` (case-insensitive).
- Plot titles in `go_enrichment.R` mention NRF1 targets; edit the script if you reuse it for another gene list.
- Large downloads and Docker image pulls need network access and disk space.
