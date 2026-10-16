
# ============================================================
# 03_functional_enrichment.R
# GO and KEGG enrichment analysis
# ============================================================

# Set the project working directory
setwd("C:/Users/ACER/OneDrive/rnaseq-project")

# Create output directories if they do not exist
dir.create("results/GO", recursive = TRUE, showWarnings = FALSE)
dir.create("results/KEGG", recursive = TRUE, showWarnings = FALSE)
dir.create("Data/Plots", recursive = TRUE, showWarnings = FALSE)

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

# Load annotated differential-expression results
annotated_results <- read.csv(
  "results/deseq2_annotated_genes.csv",
  stringsAsFactors = FALSE
)

# Inspect the imported data
cat("Rows loaded:", nrow(annotated_results), "\n")
cat("Column names:\n")
print(names(annotated_results))

# Check that essential columns are present
required_columns <- c("gene_id", "log2FoldChange", "padj")

if (!all(required_columns %in% names(annotated_results))) {
  stop("Required columns are missing. Check the column names above.")
}

cat("\nInput file loaded successfully.\n")

# ============================================================
# 1. Prepare significant, upregulated and downregulated genes
# ============================================================

# Use adjusted p-value < 0.05
sig_annotated <- subset(
  annotated_results,
  !is.na(padj) & padj < 0.05
)

# Separate upregulated and downregulated genes
up_annotated <- subset(
  sig_annotated,
  log2FoldChange > 0
)

down_annotated <- subset(
  sig_annotated,
  log2FoldChange < 0
)

cat("\nDifferential expression summary:\n")
cat("Significant genes:", nrow(sig_annotated), "\n")
cat("Upregulated genes:", nrow(up_annotated), "\n")
cat("Downregulated genes:", nrow(down_annotated), "\n")


# ============================================================
# 2. Prepare Entrez IDs for GO analysis
# ============================================================

# Convert IDs to character format, remove missing values,
# and remove duplicates
sig_entrez_ids <- unique(
  as.character(na.omit(sig_annotated$ENTREZID))
)

up_ids <- unique(
  as.character(na.omit(up_annotated$ENTREZID))
)

down_ids <- unique(
  as.character(na.omit(down_annotated$ENTREZID))
)

background_ids <- unique(
  as.character(na.omit(annotated_results$ENTREZID))
)

cat("\nGenes available for GO analysis:\n")
cat("Significant:", length(sig_entrez_ids), "\n")
cat("Upregulated:", length(up_ids), "\n")
cat("Downregulated:", length(down_ids), "\n")
cat("Background:", length(background_ids), "\n")

# Check that the IDs are character vectors
stopifnot(
  is.character(sig_entrez_ids),
  is.character(up_ids),
  is.character(down_ids),
  is.character(background_ids)
)

# ============================================================
# 3. GO Biological Process enrichment
# ============================================================

ego_bp <- enrichGO(
  gene = sig_entrez_ids,
  universe = background_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_up <- enrichGO(
  gene = up_ids,
  universe = background_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_down <- enrichGO(
  gene = down_ids,
  universe = background_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ============================================================
# 4. Export GO enrichment tables
# ============================================================

write.csv(
  as.data.frame(ego_bp),
  "results/GO/GO_BP_enrichment.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_up),
  "results/GO/GO_BP_upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_down),
  "results/GO/GO_BP_downregulated.csv",
  row.names = FALSE
)

cat("\nGO enrichment results saved.\n")
cat("Combined GO terms:", nrow(as.data.frame(ego_bp)), "\n")
cat("Upregulated GO terms:", nrow(as.data.frame(ego_up)), "\n")
cat("Downregulated GO terms:", nrow(as.data.frame(ego_down)), "\n")


# ============================================================
# 5. KEGG pathway enrichment analysis
# ============================================================

# Run KEGG enrichment for upregulated genes
kegg_up <- enrichKEGG(
  gene = up_ids,
  universe = background_ids,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Run KEGG enrichment for downregulated genes
kegg_down <- enrichKEGG(
  gene = down_ids,
  universe = background_ids,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Convert results into data frames
kegg_up_df <- as.data.frame(kegg_up)
kegg_down_df <- as.data.frame(kegg_down)

# Save results
write.csv(
  kegg_up_df,
  "results/KEGG/KEGG_upregulated.csv",
  row.names = FALSE
)

write.csv(
  kegg_down_df,
  "results/KEGG/KEGG_downregulated.csv",
  row.names = FALSE
)

# Print summary
cat("\nKEGG enrichment results saved.\n")
cat("Upregulated KEGG pathways:", nrow(kegg_up_df), "\n")
cat("Downregulated KEGG pathways:", nrow(kegg_down_df), "\n")


# ============================================================
# 6. Generate GO and KEGG enrichment plots
# ============================================================

# GO: combined significant genes
if (nrow(as.data.frame(ego_bp)) > 0) {
  p_go <- dotplot(
    ego_bp,
    showCategory = 15,
    title = "GO Biological Process Enrichment"
  )

  ggsave(
    "Data/Plots/GO_BP_dotplot.png",
    plot = p_go,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# GO: upregulated genes
if (nrow(as.data.frame(ego_up)) > 0) {
  p_go_up <- dotplot(
    ego_up,
    showCategory = 15,
    title = "GO BP Enrichment - Upregulated Genes"
  )

  ggsave(
    "Data/Plots/GO_BP_upregulated.png",
    plot = p_go_up,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# GO: downregulated genes
if (nrow(as.data.frame(ego_down)) > 0) {
  p_go_down <- dotplot(
    ego_down,
    showCategory = 15,
    title = "GO BP Enrichment - Downregulated Genes"
  )

  ggsave(
    "Data/Plots/GO_BP_downregulated.png",
    plot = p_go_down,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# KEGG: upregulated genes
if (nrow(kegg_up_df) > 0) {
  p_kegg_up <- dotplot(
    kegg_up,
    showCategory = 15,
    title = "KEGG Enrichment - Upregulated Genes"
  )

  ggsave(
    "Data/Plots/KEGG_upregulated.png",
    plot = p_kegg_up,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# KEGG: downregulated genes
if (nrow(kegg_down_df) > 0) {
  p_kegg_down <- dotplot(
    kegg_down,
    showCategory = 15,
    title = "KEGG Enrichment - Downregulated Genes"
  )

  ggsave(
    "Data/Plots/KEGG_downregulated.png",
    plot = p_kegg_down,
    width = 10,
    height = 8,
    dpi = 300
  )
}

cat("\nEnrichment plots saved.\n")
