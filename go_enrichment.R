suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else "output/nrf1_genes.txt"
out_dir <- if (length(args) >= 2) args[[2]] else "output"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genes <- readLines(input_file, warn = FALSE)
genes <- genes[genes != "" & !is.na(genes)]
cat("Total genes loaded:", length(genes), "\n")

if (length(genes) == 0) {
  stop("No genes found in input file.")
}

gene_ids <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
cat("Successfully mapped to Entrez IDs:", nrow(gene_ids), "\n")

if (nrow(gene_ids) == 0) {
  stop("No gene symbols were mapped to Entrez IDs.")
}

ego <- enrichGO(
  gene = unique(gene_ids$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

results_df <- as.data.frame(ego)
cat("Enriched GO terms found:", nrow(results_df), "\n")

write.csv(results_df, file.path(out_dir, "go_results.csv"), row.names = FALSE)
cat("Saved:", file.path(out_dir, "go_results.csv"), "\n")

if (nrow(results_df) > 0) {
  png(file.path(out_dir, "go_dotplot.png"), width = 1400, height = 900, res = 130)
  print(dotplot(ego, showCategory = 20, title = "GO BP — NRF1 target genes"))
  dev.off()
  cat("Saved:", file.path(out_dir, "go_dotplot.png"), "\n")

  png(file.path(out_dir, "go_barplot.png"), width = 1400, height = 900, res = 130)
  print(barplot(ego, showCategory = 20, title = "GO enrichment — NRF1 target genes"))
  dev.off()
  cat("Saved:", file.path(out_dir, "go_barplot.png"), "\n")

  png(file.path(out_dir, "go_cnetplot.png"), width = 1600, height = 1200, res = 130)
  print(cnetplot(ego, showCategory = min(10, nrow(results_df))))
  dev.off()
  cat("Saved:", file.path(out_dir, "go_cnetplot.png"), "\n")
}

cat("All done.\n")
