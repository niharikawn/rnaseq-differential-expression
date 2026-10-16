
# RNA-seq Differential Expression Analysis

## Project Overview

This project analyzes RNA-seq data to identify genes differentially expressed between dexamethasone-treated and untreated human airway smooth muscle cells.

The analysis uses the Airway dataset available through Bioconductor and the DESeq2 package in R.

## Objectives

- Inspect RNA-seq count data and sample metadata.
- Perform differential expression analysis.
- Visualize sample relationships using PCA.
- Generate volcano and MA plots.
- Annotate genes and identify enriched biological processes and pathways using Gene Ontology (GO) and KEGG.
- Summarize the top upregulated and downregulated genes.

## Dataset

- Dataset: Airway RNA-seq
- Organism: Homo sapiens
- Samples: 8
- Experimental groups: 4 untreated and 4 dexamethasone-treated samples
- Design formula: `~ cell + dex`

The model accounts for cell-line differences while estimating the treatment effect.

## Tools and Packages

- R
- Bioconductor
- DESeq2
- airway
- AnnotationDbi
- org.Hs.eg.db
- clusterProfiler
- enrichplot
- ggplot2

## Analysis Workflow

1. Inspect the dataset and sample metadata.
2. Filter genes with very low counts.
3. Perform differential expression analysis using DESeq2.
4. Generate PCA, volcano, and MA plots.
5. Annotate differentially expressed genes.
6. Perform GO Biological Process and KEGG enrichment analyses.
7. Export result tables and visualizations.

## Key Results

Using an adjusted p-value threshold of 0.05, the analysis identified 4,100 statistically significant genes in the DESeq2 results.

The GO and KEGG analyses identified enriched biological processes and pathways among upregulated and downregulated genes.

Enrichment results indicate statistical overrepresentation; they do not independently establish pathway activation or inhibition.


## Project Structure

```text
rnaseq-project/
├── R/
│   ├── 01_inspect_data.R
│   ├── 02_differential_expression.R
│   └── 03_functional_enrichment.R
├── Data/
│   └── Plots/
├── results/
│   ├── GO/
│   ├── KEGG/
│   ├── QC/
│   ├── dds.rds
│   ├── normalized_counts.csv
│   ├── deseq2_annotated_genes.csv
│   └── top40_differentially_expressed_genes.csv
├── README.md
├── .gitignore
└── test_setup.R
```

## Reproducibility

Run the R scripts in numerical order from the project root. Install the required R and Bioconductor packages before running the pipeline.

The differential expression script saves the fitted DESeq2 object as `results/dds.rds`.

## Limitations

The analysis includes eight samples from four cell lines. Results should be interpreted in the context of this experimental design and the limitations of the dataset.

## Author

Niharika W N
MSc Bioinformatics | Computational Biology | Genomics

GitHub: @niharikawn

LinkedIn: Niharika Wooluvarana Naveen

ORCID: 0009-0004-6345-7064
