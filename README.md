## ⚠️ This project is currently under development. Some features may be incomplete or subject to change.
# EdgeR Analysis for Colon Cancer Gene Expression

## Overview
This project implements differential gene expression analysis using the `edgeR` package in R, focusing on gene expression data related to colon cancer. The analysis compares gene expression between patients based on their survival outcomes, aiming to identify significant genes and pathways involved in colon cancer.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Analysis Steps](#analysis-steps)
- [Results](#results)
 
  

## Prerequisites
Ensure you have R installed along with the necessary packages. This project uses:
- `tidyverse`
- `edgeR`
- `fgsea`
- `msigdbr`
- `Cairo`

## Installation
To install the required R packages, run the following commands in your R console:

```R
install.packages(c("tidyverse", "edgeR", "fgsea", "msigdbr", "Cairo"))
```
## Analysis Steps

1. **Load Necessary Libraries**: Load all required libraries.
2. **Load Data**: Load the colon cancer dataset (`ex3.Rdata`).
3. **Data Preparation**: Filter the counts and survival data for overlapping samples.
4. **Create DGEList Object**: Convert the counts data to a DGEList object and filter lowly expressed genes.
5. **Normalization**: Normalize the count data.
6. **Design Matrix Setup**: Create a design matrix for the survival analysis.
7. **Model Fitting**: Estimate dispersion and fit the model using `glmQLFit`.
8. **Differential Expression Analysis**: Conduct the analysis using `glmQLFTest`.
9. **Visualizations**: Create a volcano plot to visualize the results.
10. **Gene Ranking**: Prepare a ranked list of genes based on log fold change.
11. **GSEA Analysis**: Load gene sets and perform GSEA using the ranked gene list.
12. **Top Pathways Visualization**: Create a bar plot for the top enriched pathways.
13. **Export Results**: Save the expression list to CSV and TXT files.

## Results

The results include a volcano plot visualizing the differential expression analysis and a bar plot showing the top enriched pathways. Additionally, a CSV file containing the significant genes and their corresponding log fold changes and adjusted p-values is generated.

