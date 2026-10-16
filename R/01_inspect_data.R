
# RNA-seq Differential Expression Pipeline
# Script 01: Inspect the example dataset

# Load the required package
library(airway)
library(SummarizedExperiment)

# Load the example RNA-seq dataset
data("airway")

# Extract the gene-level count matrix
counts <- assay(airway)

# Extract sample metadata
sample_metadata <- as.data.frame(colData(airway))

# Print basic dataset information
cat("Number of genes/features:", nrow(counts), "\n")
cat("Number of samples:", ncol(counts), "\n")

cat("\nSample names:\n")
print(colnames(counts))

cat("\nCount matrix preview:\n")
print(counts[1:6, 1:4])

cat("\nSample metadata:\n")
print(sample_metadata)

# Verify that counts and metadata match
stopifnot(ncol(counts) == nrow(sample_metadata))

cat("\nSUCCESS: Count matrix and sample metadata match.\n")
