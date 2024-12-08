# **Biomarker and Signaling Pathway Discovery to Improve Patient Survival in Colon Cancer**

## **Overview**
This repository contains code and data analysis pipelines for discovering biomarkers and signaling pathways to improve the survival of patients with colon cancer. The project leverages RNA-seq data to perform differential expression analysis, gene set enrichment analysis (GSEA), single-sample GSEA (ssGSEA), and survival analysis to identify and validate potential biomarkers.

## **Features**
- **RNA-seq Differential Expression Analysis** using the `edgeR` package.
- **Gene Set Enrichment Analysis (GSEA)** for identifying enriched pathways.
- **Single-Sample GSEA (ssGSEA)** to analyze individual pathway activity scores.
- **Kaplan-Meier Survival Analysis** to evaluate survival differences.
- Visualizations including:
  - Volcano plots for significant genes.
  - Bar plots for top enriched pathways.
  - Kaplan-Meier survival curves.

---

## **Requirements**
This project uses R and the following packages:
- `tidyverse`
- `edgeR`
- `limma`
- `ggplot2`
- `ggrepel`
- `fgsea`
- `msigdbr`
- `Cairo`
- `GSVA`
- `survival`

Ensure you have R installed on your system. You can install the required packages using the following command in R:
```R
install.packages(c("tidyverse", "edgeR", "limma", "ggplot2", "ggrepel", "fgsea", "msigdbr", "Cairo", "GSVA", "survival"))
