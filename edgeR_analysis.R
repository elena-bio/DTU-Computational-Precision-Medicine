## Load necessary packages
library(tidyverse)
library(edgeR)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(Cairo)
library(ggrepel)

# Look day 7: for Kaplan Meier plots
# Look day 4

# Load the data
load("/home/projects/22123/ex3.Rdata")

# -------------------------------------------
# Data Preparation
# -------------------------------------------

# Step 1: Filter the counts data and the survival data for overlapping samples
counts_samples <- unique(COAD_counts$sample)
surv_samples <- unique(surv$sample)
overlapping_samples <- intersect(counts_samples, surv_samples)

COAD_counts_filtered <- filter(COAD_counts, sample %in% overlapping_samples)
surv_filtered <- filter(surv, sample %in% overlapping_samples)

# Step 2: Convert the counts data to a matrix with genes as rows and samples as columns
COAD_counts_matrix <- COAD_counts_filtered |> 
  pivot_wider(names_from = sample, values_from = count) |> 
  column_to_rownames("gene")

# Step 3: Prepare the survival data for analysis
surv_filtered <- column_to_rownames(surv_filtered, "sample")

# Step 4: Ensure column names in the counts matrix match row names in the survival data
COAD_counts_matrix <- COAD_counts_matrix[, overlapping_samples]
surv_filtered <- surv_filtered[overlapping_samples,]

# Step 5: Categorize patients into long-term and short-term survivors based on median OS.time
median_OS_time <- median(surv_filtered$OS.time)
surv_filtered$group <- ifelse(surv_filtered$OS.time <= median_OS_time, "Short-term", "Long-term")

# -------------------------------------------
# EdgeR Differential Expression Analysis
# -------------------------------------------

# Step 6: Create a DGEList object and filter lowly expressed genes
dge <- DGEList(counts = COAD_counts_matrix)
keep <- filterByExpr(dge)
dge <- dge[keep, keep.lib.sizes = FALSE]

# Step 7: Normalize the data
dge <- calcNormFactors(dge)

# Step 8: Set up the design matrix for survival analysis
design <- model.matrix(~ surv_filtered$OS)
colnames(design) <- c("Intercept", "OS")

# Step 9: Estimate dispersion and fit the model
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Step 10: Conduct differential expression analysis
qlf <- glmQLFTest(fit, coef = "OS")
edgeR_res <- topTags(qlf, n = Inf)

# Step 11: Convert results to a data frame and identify significant genes
edgeR_res_df <- as.data.frame(edgeR_res) |> 
  mutate(significant = FDR < 0.05)

# Step 12: Select top upregulated and downregulated genes for labeling
top_genes <- edgeR_res_df %>%
  dplyr::filter(significant == TRUE) %>%
  dplyr::arrange(desc(logFC)) %>%
  dplyr::slice(1:10)

bottom_genes <- edgeR_res_df %>%
  dplyr::filter(significant == TRUE) %>%
  dplyr::arrange(logFC) %>%
  dplyr::slice(1:10)

extreme_genes <- bind_rows(top_genes, bottom_genes)
extreme_genes <- extreme_genes %>% rownames_to_column(var = "gene_names")

# Step 13: Create a volcano plot
ggplot(edgeR_res_df, aes(x = logFC, y = -log10(FDR), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot (edgeR)") + 
  theme_minimal() +
  geom_text_repel(data = extreme_genes, aes(label = gene_names), 
                  max.overlaps = 20, box.padding = 0.2, 
                  point.padding = 0.1, segment.color = "grey50")

# Step 14: Prepare ranked list of genes for GSEA
ranked_genes <- edgeR_res_df$logFC
names(ranked_genes) <- rownames(edgeR_res_df)
ranked_genes <- ranked_genes[!is.na(ranked_genes)]

# -------------------------------------------
# GSEA Analysis
# -------------------------------------------

# Step 15: Load gene sets for fgsea
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)

# Step 16: Run fgsea
fgsea_res <- fgsea(pathways = gene_sets, stats = ranked_genes, minSize = 15, maxSize = 500)

# Step 17: Create a bar plot showing the top 5 enriched pathways
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

# -------------------------------------------
# Kaplan-Meier Analysis
# -------------------------------------------

# Step 18: Kaplan-Meier survival curves for short-term vs. long-term survivors

# -------------------------------------------
# Save Results
# -------------------------------------------

# Step 19: Save differential expression results

# Step 20: Save GSEA results

