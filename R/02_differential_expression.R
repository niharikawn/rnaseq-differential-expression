
# ============================================================
# RNA-seq Differential Expression Pipeline
# Script 02: Differential Expression with DESeq2
# Dataset: Airway RNA-seq (dexamethasone treatment)
# ============================================================

# Step 1: Load required packages
library(airway)
library(DESeq2)
library(ggplot2)

# Step 2: Load the example dataset
data("airway")

# Extract count matrix and sample metadata
counts_matrix <- assay(airway)
sample_metadata <- as.data.frame(colData(airway))

# Set untreated samples as the reference group
sample_metadata$dex <- relevel(
  factor(sample_metadata$dex),
  ref = "untrt"
)

# Account for cell-line differences and treatment
design_formula <- ~ cell + dex

# Create output directories
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("Data/Plots", showWarnings = FALSE, recursive = TRUE)

# Step 3: Build the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),
  colData = sample_metadata,
  design = design_formula
)

# Filter genes with very low counts
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

cat("Genes retained after filtering:", nrow(dds), "\n")

# Step 4: Run differential-expression analysis
dds <- DESeq(dds)

# Save the fitted DESeq2 object for future use
saveRDS(dds, file = "results/dds.rds")

cat("DESeq2 object saved to results/dds.rds\n")

# Step 5: Compare treated samples against untreated samples
res <- results(
  dds,
  contrast = c("dex", "trt", "untrt"),
  alpha = 0.05
)

# Sort results by adjusted p-value
res_ordered <- res[order(res$padj, na.last = TRUE), ]

# Display the top results
print(head(as.data.frame(res_ordered), 10))

# Save the complete differential-expression results
results_df <- as.data.frame(res_ordered)
results_df$gene_id <- rownames(results_df)

write.csv(
  results_df,
  file = "results/differential_expression_results.csv",
  row.names = FALSE
)

cat("Differential-expression results saved.\n")

# Step 6: Create a volcano plot
volcano_data <- as.data.frame(res)
volcano_data$gene_id <- rownames(volcano_data)

# Classify genes using adjusted p-value and fold-change thresholds
volcano_data$status <- "Not significant"

volcano_data$status[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange >= 1
] <- "Upregulated"

volcano_data$status[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange <= -1
] <- "Downregulated"

# Prevent infinite values when adjusted p-values are zero
volcano_data$padj_plot <- pmax(volcano_data$padj, 1e-300)

p_volcano <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj_plot),
    color = status
  )
) +
  geom_point(alpha = 0.65, size = 1.5, na.rm = TRUE) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "RNA-seq Differential Expression",
    subtitle = "Dexamethasone-treated vs untreated",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value",
    color = "Expression status"
  ) +
  theme_minimal(base_size = 12)

print(p_volcano)

ggsave(
  filename = "Data/Plots/volcano_plot.png",
  plot = p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Volcano plot saved.\n")

# Step 7: Create an MA plot
png(
  filename = "Data/Plots/MA_plot.png",
  width = 2400,
  height = 1800,
  res = 300
)

plotMA(
  res,
  alpha = 0.05,
  main = "MA Plot: Treated vs Untreated",
  ylim = c(-5, 5)
)

dev.off()

cat("MA plot saved.\n")

# Step 8: Principal Component Analysis (PCA)
# Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Extract PCA data using treatment and cell-line metadata
pca_data <- plotPCA(
  vsd,
  intgroup = c("dex", "cell"),
  returnData = TRUE
)

# Calculate percentage of variance explained
percent_var <- round(
  100 * attr(pca_data, "percentVar"),
  1
)

p_pca <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = dex,
    shape = cell
  )
) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of RNA-seq Samples") +
  theme_minimal(base_size = 12)

print(p_pca)

ggsave(
  filename = "Data/Plots/PCA_plot.png",
  plot = p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

cat("PCA plot saved.\n")

# Step 9: Report analysis summary
cat("\n====================================\n")
cat("Differential-expression analysis complete\n")
cat("Design: ~ cell + dex\n")
cat("Contrast: dex trt vs untrt\n")
cat("Genes tested:", sum(!is.na(res$pvalue)), "\n")
cat(
  "Significant genes (adjusted p < 0.05):",
  sum(res$padj < 0.05, na.rm = TRUE),
  "\n"
)
cat("Saved DESeq2 object: results/dds.rds\n")
cat("====================================\n")
