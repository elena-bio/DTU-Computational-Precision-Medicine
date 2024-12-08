## Load necessary packages
library(tidyverse)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(Cairo)
library(GSVA)
library(survival)  #For Kaplan-Meier analysis

## NOTE: This script covers the key steps for analyzing RNA-seq .
## including differential expression analysis, GSEA, ssGSEA, and survival analysis.


## Load the data for this week's exercises
load("/home/projects/22123/ex3.Rdata")

# -------------- Data Preparation -------------- #

# Step 1: Filter the counts data and the survival data for overlapping(common) samples
counts_samples <- unique(COAD_counts$sample)

#surv_samples <- surv$sample
surv_samples <- unique(surv$sample) 

# Step 2: Find the common data of both objects
overlapping_samples <- intersect(counts_samples, surv_samples)

# Step 3: Filter the counts and survival data to include just common samples  
COAD_counts_filtered <- filter(COAD_counts, sample %in% overlapping_samples)

surv_filtered <- filter(surv, sample %in% overlapping_samples)

# Step 4: Convert the counts data to a matrix with genes as rows and samples as columns
COAD_counts_matrix <- COAD_counts_filtered |> 
  pivot_wider(names_from = sample, values_from = count) |> 
  column_to_rownames("gene")

# Step 5: # Convert the "sample" column to row names for compatibility with analysis
surv_filtered <- column_to_rownames(surv_filtered, "sample")

# Step6: Make sure column names in the counts matrix match row names in the survival data
COAD_counts_matrix <- COAD_counts_matrix[, overlapping_samples]
surv_filtered <- surv_filtered[overlapping_samples,]

# Step 7: Categorize patients into long-term and short-term survivors based on median OS.time
median_OS_time <- median(surv_filtered$OS.time)

# Step 8:Split patients into two groups: short-term and long-term survivors based on median survival time
surv_filtered$group <- ifelse(surv_filtered$OS.time > median_OS_time, "Long-term", "Short-term")

# View the distribution of patients in each group
table(surv_filtered$group)

# -------------- Differential Expression Analysis with edgeR -------------- #

# Step9: Create a DGEList object and filter lowly expressed genes
dge <- DGEList(counts = COAD_counts_matrix, group = surv_filtered$group)

#Step 10: Filter out lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, keep.lib.sizes = FALSE]

# Step 11: Normalize the data
dge <- calcNormFactors(dge, method = "TMM")

#Check normalization
#head(dge$samples, 2)

# Step 12: Set up the design matrix for survival analysis with edgeR
design <- model.matrix(~ surv_filtered$group)
colnames(design) <- c("Intercept", "Group")

# Step 13: Estimate dispersion and fit model
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Step 14: Conduct differential expression analysis
qlf <- glmQLFTest(fit, coef = "Group")
edgeR_res <- topTags(qlf, n = Inf)

#head(edgeR_res$table)

# Convert the results to a data frame and add significance column based on adjusted p-value
edgeR_res_df <- as.data.frame(edgeR_res) |> 
  mutate(significant = FDR < 0.05)

# Step 15: Filter significant genes
significant_genes <- edgeR_res_df[edgeR_res_df$FDR < 0.05, ]

# -------------- Visualization: Volcano Plot -------------- #

ggplot(edgeR_res_df, 
       aes(x = logFC, 
           y = -log10(FDR), 
           color = logFC, 
           label = rownames(edgeR_res_df))) + 
  geom_point() +  
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +  
  geom_text_repel(data = significant_genes,  # Use only significant genes for labels
                  aes(label = rownames(significant_genes)), 
                  size = 3, 
                  box.padding = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey') +  
  labs(x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value", 
       title = "Volcano Plot (Significant Genes)") + 
  theme_minimal() +
  theme(legend.position = "right")

# -------------- Pathway Analysis: GSEA -------------- #

# Step 16:(GSEA) Prepare ranked list of genes by test statistic
ranked_genes <- edgeR_res_df$logFC
names(ranked_genes) <- rownames(edgeR_res_df)

# Step 17:Remove missing values from the ranked gene list
ranked_genes <- ranked_genes[!is.na(ranked_genes)]

# Step 18: Load gene sets for fgsea analysis
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)

# Step 19: Run fgsea using the gene sets and the ranked list of genes
fgsea_res <- fgsea(pathways = gene_sets, stats = ranked_genes, minSize = 15, maxSize = 500)
#head(fgsea_res, 2)

# Step 20: Filter for pathways with positive NES (upregulated)
topUpregulatedPathways <- fgsea_res %>%
  dplyr::filter(NES > 0) %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::slice(1:5)

# Create the bar plot for the top 5 upregulated pathways
ggplot(topUpregulatedPathways, 
       mapping = aes(x = reorder(pathway, NES), 
                     y = NES)) +
  geom_col() +
  coord_flip() +  # Flip the axes
  labs(x = "Pathway", 
       y = "Normalized Enrichment Score (NES)", 
       title = "Top 5 Upregulated Pathways") +
  theme_minimal()

# Step 21: Output the expression list with additional information
expression_list <- edgeR_res_df %>%
  mutate(gene = rownames(edgeR_res_df)) %>%
  select(gene, log2FoldChange = logFC, padj = FDR, significant)


# -------------- ssGSEA Analysis --------------

# Step 1: Normalize counts to TPM
# Ensure counts are converted to TPM format
counts_tpm <- apply(COAD_counts_matrix, 2, function(x) {
  if (sum(x) > 0) {  # Prevent division by zero
    (x / sum(x)) * 1e6
  } else {
    rep(0, length(x))
  }
})

# Step 2: Load gene sets for ssGSEA analysis
library(msigdbr)
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)

# Step 3: Identify top pathways using fgsea
# Ensure ranked_genes is provided as a named numeric vector
library(fgsea)
if (!is.numeric(ranked_genes) || is.null(names(ranked_genes))) {
  stop("Error: `ranked_genes` must be a named numeric vector.")
}

# Run fgsea to rank pathways
fgsea_res <- fgsea(pathways = gene_sets, stats = ranked_genes, minSize = 15, maxSize = 500)

# Select the top 10 most enriched pathways based on NES
library(dplyr)
topPathways <- fgsea_res %>%
  arrange(desc(abs(NES))) %>%
  slice(1:10)  # Select top 10 pathways

# Extract the names of the top pathways
selected_pathways <- topPathways$pathway

# Step 4: Create the list of gene sets based on selected pathways
selected_gene_sets <- gene_sets[selected_pathways]

# Step 5: Run ssGSEA using TPM and selected pathways
library(GSVA)
ssgsea_results <- gsva(
  expr = counts_tpm,             # TPM-normalized counts
  gset.idx.list = selected_gene_sets,  # Top enriched pathways
  method = "ssgsea",             # ssGSEA method
  parallel.sz = 1                # Number of parallel workers
)

# Step 6: Convert ssGSEA results to a data frame
ssgsea_results_df <- as.data.frame(t(ssgsea_results))
ssgsea_results_df$sample <- rownames(ssgsea_results_df)

# Step 7: Visualize Enrichment Score Distributions
library(tidyr)
library(ggplot2)

# Transform ssGSEA results to long format for ggplot
long_data <- ssgsea_results_df %>%
  pivot_longer(
    cols = all_of(selected_pathways),  # Use selected pathways
    names_to = "pathway",
    values_to = "EnrichmentScore"
  )

# Create the histogram plot
combined_histogram_plot <- ggplot(long_data, aes(x = EnrichmentScore)) +
  geom_histogram(binwidth = 0.02, fill = "lightblue", color = "black") +
  facet_wrap(~ pathway, scales = "free_y", ncol = 3) +  # Grid layout for pathways
  labs(
    title = "Top 10 Enriched Pathways - ssGSEA",  # Updated title
    x = "Enrichment Score",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    strip.text = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# Display the plot
print(combined_histogram_plot)

# -------------- survival Analysis --------------

# Ensure the survival dataset `surv_filtered` has a `sample` column
if (!"sample" %in% colnames(surv_filtered)) {
  library(tibble)
  surv_filtered <- surv_filtered %>%
    rownames_to_column(var = "sample")  # Convert rownames to a `sample` column
}

# Ensure the ssGSEA results `ssgsea_results_df` has a `sample` column
if (!"sample" %in% colnames(ssgsea_results_df)) {
  ssgsea_results_df$sample <- rownames(ssgsea_results_df)  # Add rownames as `sample` column
}

# Align sample names (remove any extra spaces)
surv_filtered$sample <- trimws(surv_filtered$sample)
ssgsea_results_df$sample <- trimws(ssgsea_results_df$sample)

# Filter data frames to include only common samples
common_samples <- intersect(surv_filtered$sample, ssgsea_results_df$sample)
surv_filtered <- surv_filtered %>% filter(sample %in% common_samples)
ssgsea_results_df <- ssgsea_results_df %>% filter(sample %in% common_samples)

# Rename OS to OS.event for consistency if it exists
if ("OS" %in% colnames(surv_filtered)) {
  surv_filtered <- surv_filtered %>%
    rename(OS.event = OS)
}

# Merge survival data with ssGSEA results
merged_data <- left_join(surv_filtered, ssgsea_results_df, by = "sample")

# Ensure no missing OS.time or OS.event in merged_data
merged_data <- merged_data %>%
  filter(!is.na(OS.time) & !is.na(OS.event))

# Check the final dataset structure
#cat("Head of merged_data:\n")
#print(head(merged_data))
#cat("\nDimensions of merged_data:\n")
#print(dim(merged_data))

# Remove rows with missing OS.time or OS.event
merged_data <- merged_data %>%
  filter(!is.na(OS.time) & !is.na(OS.event))

# Split patients based on the enrichment score of the first selected pathway
merged_data$group <- ifelse(
  merged_data[["GOBP_VENTRICULAR_CARDIAC_MUSCLE_TISSUE_MORPHOGENESIS"]] > median(merged_data[["GOBP_VENTRICULAR_CARDIAC_MUSCLE_TISSUE_MORPHOGENESIS"]]),
  "High",
  "Low"
)
# Create a Surv object
#library(survival)
# surv_object <- Surv(time = merged_data$OS.time, event = merged_data$OS.event)

# Fit the Kaplan-Meier model
# km_fit <- survfit(surv_object ~ group, data = merged_data)

# Create the Kaplan-Meier plot
library(ggplot2)
library(survival)
km_fit <- survfit(Surv(time = merged_data$OS.time, event = merged_data$OS.event) ~ group, data = merged_data)
#str(km_fit)

# Create a data frame from the Kaplan-Meier fit
km_data <- data.frame(
  time = km_fit$time,         # Survival times
  surv = km_fit$surv,         # Survival probabilities
  group = rep(names(km_fit$strata), times = km_fit$strata),  # Group labels
  n.risk = km_fit$n.risk,     # Number of patients at risk
  n.event = km_fit$n.event    # Number of events
)

#head(km_data)

library(ggplot2)

km_plot <- ggplot(data = km_data, aes(x = time, y = surv, color = group)) +
  geom_step(size = 1) +  # Step function for survival probabilities
  labs(
    title = "Kaplan-Meier Survival Curve",
    x = "Time (Months)",
    y = "Survival Probability",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Display the Kaplan-Meier plot
print(km_plot)

# Perform a log-rank test to evaluate survival differences between groups
logrank_test <- survdiff(Surv(time = merged_data$OS.time, event = merged_data$OS.event) ~ group, data = merged_data)

# View the results of the log-rank 
print(logrank_test)















































































































































































