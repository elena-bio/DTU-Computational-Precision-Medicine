##edgeR analysis
## Load necessary packages 
library(tidyverse)
library(edgeR)
library(fgsea)
library(msigdbr)
install.packages("Cairo")
library(Cairo)

setwd("/net/pupil1/home/people/s232130/TidayThusday/edger")

## Load the data for this week's exercises
load("ex3.Rdata")

## Step 1: Filter the counts data and the survival data for overlapping samples
counts_samples <- unique(COAD_counts$sample)
surv_samples <- surv$sample
overlapping_samples <- intersect(counts_samples, surv_samples)

COAD_counts_filtered <- filter(COAD_counts, sample %in% overlapping_samples)
surv_filtered <- filter(surv, sample %in% overlapping_samples)

## Step 2: Convert the counts data to a matrix with genes as rows and samples as columns
COAD_counts_matrix <- COAD_counts_filtered |> 
  pivot_wider(names_from = sample, values_from = count) |> 
  column_to_rownames("gene")

## Step 3: Prepare the survival data for analysis
surv_filtered <- column_to_rownames(surv_filtered, "sample")

## Make sure column names in the counts matrix match row names in the survival data
COAD_counts_matrix <- COAD_counts_matrix[, overlapping_samples]
surv_filtered <- surv_filtered[overlapping_samples,]

## Step 4: Create a DGEList object and filter lowly expressed genes
dge <- DGEList(counts = COAD_counts_matrix)

# Filter out lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, keep.lib.sizes = FALSE]

## Step 5: Normalize the data
## Step 5: Normalize the data
dge <- calcNormFactors(dge)

## Step 6: Set up the design matrix for survival analysis with edgeR
design <- model.matrix(~ surv_filtered$OS)
colnames(design) <- c("Intercept", "OS")

## Step 7: Estimate dispersion and fit model
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

## Step 8: Conduct differential expression analysis
qlf <- glmQLFTest(fit, coef = "OS")
edgeR_res <- topTags(qlf, n = Inf)

## Convert the results to a data frame and add significance column based on adjusted p-value
edgeR_res_df <- as.data.frame(edgeR_res) |> 
  mutate(significant = FDR < 0.05)

## Step 9: Create a volcano plot
ggplot(edgeR_res_df, 
       aes(x = logFC, 
           y = -log10(FDR), 
           color = significant)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value", 
       title = "Volcano Plot (edgeR)") + 
  theme_minimal()

## Step 10: Prepare ranked list of genes by test statistic
ranked_genes <- edgeR_res_df$logFC
names(ranked_genes) <- rownames(edgeR_res_df)

## Remove missing values from the ranked gene list
ranked_genes <- ranked_genes[!is.na(ranked_genes)]

## Step 11: Load gene sets for fgsea analysis
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)

## Step 12: Run fgsea using the gene sets and the ranked list of genes
## Step 12: Run fgsea using the gene sets and the ranked list of genes
# Ensure ranked_genes is named properly
ranked_genes <- edgeR_res_df$logFC  # Use logFC from edgeR results for ranking
names(ranked_genes) <- rownames(edgeR_res_df)  # Assign gene names as names

# Remove missing values from the ranked gene list
ranked_genes <- ranked_genes[!is.na(ranked_genes)]

# Check that ranked_genes is named
if (any(is.na(names(ranked_genes)))) {
  stop("The ranked_genes vector must be named with gene identifiers.")
}

# Load gene sets for fgsea analysis
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)

# Run fgsea using the gene sets and the ranked list of genes
fgsea_res <- fgsea(pathways = gene_sets, stats = ranked_genes, minSize = 15, maxSize = 500)



## Step 13: Create a bar plot showing the normalized enrichment scores of the top 5 most enriched significant pathways 
topPathways <- fgsea_res |> 
  dplyr::arrange(padj) |> 
  dplyr::slice(1:5)

ggplot(topPathways, 
       mapping = aes(x = reorder(pathway, NES), 
                     y = NES)) +
  geom_col() +
  coord_flip() +
  labs(x = "Pathway", 
       y = "Normalized Enrichment Score (NES)", 
       title = "Top 5 Enriched Pathways") +
  theme_minimal()

## Step 14: Output the expression list with additional information
## Step 14: Output the expression list with additional information
# Create a new data frame with gene names as a column
# Create the expression list with appropriate columns
expression_list <- edgeR_res_df %>%
  mutate(gene = rownames(edgeR_res_df)) %>%  # Add a new column for gene names
  select(gene, log2FoldChange = logFC, padj = FDR, significant)

# Save the expression list to a CSV file
write.csv(expression_list, file = "expression_list_edgeR.csv", row.names = FALSE)


## Optional: Print a confirmation message
print("File saved as expression_list_edgeR.csv in the home directory.")

write.table(expression_list, file = "expression_list_edgeR.txt", row.names = FALSE, sep = "\t")

plot(expression_list$column1, expression_list$column2)  # Adjust as necessary


