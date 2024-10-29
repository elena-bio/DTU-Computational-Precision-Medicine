# EdgeR Analysis for Colon Cancer Gene Expression

## Overview
This project implements differential gene expression analysis using the `edgeR` package in R, focusing on gene expression data related to colon cancer. The analysis compares gene expression between patients based on their survival outcomes, aiming to identify significant genes and pathways involved in colon cancer.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Analysis Steps](#analysis-steps)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)

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
