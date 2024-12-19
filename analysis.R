# Install and Load Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install necessary packages
 BiocManager::install(c("GEOquery", "limma", "ggplot2", "pheatmap", "RColorBrewer"))

# Load libraries
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Download GSE10072 Dataset
cat("Downloading GSE10072 dataset...\n")
gse <- getGEO("GSE10072", GSEMatrix = TRUE)  # Download GEO dataset
data <- exprs(gse[[1]])                      # Extract expression data
meta <- pData(gse[[1]])                      # Extract metadata

# prepare Metadata and Conditions
cat("Preparing metadata...\n")
# Assign condition based on source_name_ch1
meta$condition <- factor(ifelse(meta$source_name_ch1 == "Adenocarcinoma of the Lung", "cancer", "normal"),
                         levels = c("normal", "cancer"))

# Validate the condition mapping
cat("Condition distribution:\n")
print(table(meta$condition, useNA = "ifany"))

colnames(data) <- rownames(meta)  # Align column names of data with metadata

#  Normalize Data (If Necessary)
# limma can handle log-transformed data; ensure the data is appropriate
cat("Checking data transformation...\n")
if (max(data) > 50) {  # If data is not log-transformed, apply log2 transformation
  data <- log2(data + 1)
  cat("Data log-transformed.\n")
}

# Create Design Matrix
group <- meta$condition  # Experimental groups (e.g., cancer vs normal)
design <- model.matrix(~ group)

# Fit Linear Model with limma
fit <- lmFit(data, design)  # Fit linear model
fit <- eBayes(fit)          # Empirical Bayes moderation

# Extract Differentially Expressed Genes (DEGs)
results <- topTable(fit, coef = 2, adjust = "fdr", sort.by = "P", number = Inf)  # FDR-adjusted p-values
write.csv(results, "limma_results.csv")  # Save DEGs to a CSV file

#Volcano Plot
results$logP <- -log10(results$P.Value)  # Compute -log10 p-values
results$Significant <- results$adj.P.Val < 0.05  # Mark significant genes

volcano <- ggplot(results, aes(x = logFC, y = logP)) +
  geom_point(aes(color = Significant), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-Value")
ggsave("volcano_plot_limma.png", volcano)  # Save volcano plot
print(volcano)

# Heatmap of Top DEGs
# Select top 50 DEGs by adjusted p-value
top_genes <- rownames(results[1:50, ])
heatmap_data <- data[top_genes, ]

# Scale data for better visualization
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# Create heatmap
pheatmap(heatmap_data_scaled,
         annotation_col = data.frame(Condition = group),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(brewer.pal(9, "Blues"))(50),
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename = "heatmap_limma.png")

cat("Analysis complete! Results saved as 'limma_results.csv', 'volcano_plot_limma.png', and 'heatmap_limma.png'.\n")
