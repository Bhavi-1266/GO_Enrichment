#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="${ROOT}/.venv/bin/python3"
if [[ ! -x "$PY" ]]; then PY="python3"; fi

GTF_URL="${GTF_URL:-https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz}"
HG38_URL="${HG38_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz}"

mkdir -p input output

if [[ ! -f input/gencode.v45.basic.annotation.gtf.gz ]]; then
  echo "Downloading GENCODE GTF..."
  wget -q --show-progress -O input/gencode.v45.basic.annotation.gtf.gz "$GTF_URL"
fi

if [[ ! -f input/human_gene_annotation.tsv.gz ]]; then
  echo "Building human_gene_annotation.tsv.gz from GTF..."
  "$PY" scripts/gtf_to_human_gene_annotation.py \
    --input input/gencode.v45.basic.annotation.gtf.gz \
    --output input/human_gene_annotation.tsv.gz \
    --types protein_coding
fi

echo "Building TSS BED..."
"$PY" scripts/make_tss_bed.py \
  --input input/human_gene_annotation.tsv.gz \
  --output output/genes_tss_clean.bed

if [[ ! -f input/hg38.fa.gz ]]; then
  echo "Downloading UCSC hg38 (large, one-time)..."
  wget -c --show-progress -O input/hg38.fa.gz "$HG38_URL"
fi

if [[ ! -f input/hg38.fa ]]; then
  echo "Decompressing UCSC hg38 FASTA to plain .fa for pyfaidx..."
  gzip -dc input/hg38.fa.gz > input/hg38.fa
fi

echo "Extracting promoter FASTA (creates .fai index on first run)..."
"$PY" scripts/extract_promoter_fasta.py \
  --bed output/genes_tss_clean.bed \
  --genome input/hg38.fa \
  --output output/promoters.fa

echo "Scanning NRF1 motif in promoters..."
"$PY" scripts/find_nrf1_genes.py \
  --input output/promoters.fa \
  --output output/nrf1_genes.txt

echo "Running GO enrichment in Docker (installs Bioconductor packages on first run)..."
docker run --rm \
  -v "$ROOT:/work" \
  -w /work \
  ghcr.io/bioconductor/bioconductor_docker:RELEASE_3_19 \
  bash -lc 'Rscript -e "BiocManager::install(c(\"clusterProfiler\",\"org.Hs.eg.db\"), ask=FALSE, update=FALSE, Ncpus=2)" && Rscript go_enrichment.R output/nrf1_genes.txt output'

echo "Done. See output/ for go_results.csv and plots."
