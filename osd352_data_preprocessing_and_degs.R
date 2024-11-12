#####################################
# Installing & Importing Dependencies
#####################################

if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(scales)
library(rstudioapi)

###############################
# Data Loading & Pre-processing
###############################

# Setting dynamic working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the data
star_data <- read.csv("./data/GLDS-352_rna_seq_STAR_Unnormalized_Counts.csv")

# Define conditions based on sample names
condition <- ifelse(grepl("F", colnames(star_data)[-1]), "Case", "Control")

# Preprocess function for data cleaning
preprocess_rnaseq <- function(star_data) {
  gene_id_col <- colnames(star_data)[1]
  expressed_genes <- star_data[rowSums(star_data[,-1] > 0) >= (ncol(star_data) - 1) * 0.25, ]
  log2_counts <- log2(expressed_genes[,-1] + 1)
  log2_counts[[gene_id_col]] <- expressed_genes[[gene_id_col]]
  return(list(
    raw = expressed_genes,
    log2 = log2_counts,
    gene_id_col = gene_id_col,
    filtered_genes = nrow(star_data) - nrow(expressed_genes)
  ))
}

# Apply preprocessing
processed_data <- preprocess_rnaseq(star_data)

##################################
# Differential Expression Analysis
##################################

# Prepare condition factor with proper reference level
condition <- factor(
  ifelse(grepl("F", colnames(star_data)[-1]), "Case", "Control"),
  levels = c("Control", "Case")  # Explicitly set Control as reference level
)

# Prepare DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(star_data[,-1]),
  colData = data.frame(condition = condition),
  design = ~ condition
)
rownames(dds) <- star_data[[1]]  # Set gene IDs as rownames

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results with LFC between Case and Control
res <- results(dds)  # No need to specify contrast since we set the reference level

# Order results by adjusted p-value
res_ordered <- res[order(res$padj),]

# Save differential expression results
if (!dir.exists("analysis_output")) {
  dir.create("analysis_output")
}

write.csv(as.data.frame(res_ordered), 
          "analysis_output/deseq2_results.csv", 
          row.names = TRUE)

###############################
# Visualizations
###############################

# Expression density plot
plot_expression_density <- function(log2_counts, gene_id_col) {
  melted_counts <- log2_counts %>%
    select(-all_of(gene_id_col)) %>%
    gather(key = "Sample", value = "Log2_Counts")
  
  ggplot(melted_counts, aes(x = Log2_Counts, color = Sample)) +
    geom_density(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Distribution of Log2 Expression Values",
         x = "Log2(Counts + 1)",
         y = "Density") +
    theme(legend.position = "none")
}

# Sample correlation heatmap
plot_sample_correlation <- function(log2_counts, gene_id_col) {
  cor_matrix <- cor(log2_counts %>% select(-all_of(gene_id_col)))
  pheatmap(cor_matrix,
           main = "Sample Correlation Heatmap",
           show_rownames = TRUE,
           show_colnames = TRUE,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 8,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
}

# MA Plot
plot_ma <- function(deseq_results) {
  ggplot(
    as.data.frame(deseq_results),
    aes(x = baseMean, y = log2FoldChange, color = padj < 0.05)
  ) +
    geom_point(alpha = 0.6) +
    scale_x_log10() +
    scale_color_manual(values = c("grey60", "red")) +
    labs(
      title = "MA Plot",
      x = "Mean Expression",
      y = "Log2 Fold Change",
      color = "Significant"
    ) +
    theme_minimal()
}

# Generate plots
pdf("analysis_output/quality_control_plots.pdf")
plot_expression_density(processed_data$log2, processed_data$gene_id_col)
plot_sample_correlation(processed_data$log2, processed_data$gene_id_col)
plot_ma(res_ordered)
dev.off()

###########################
# Save Processed Data
###########################

write.csv(processed_data$raw, 
          "analysis_output/processed_counts.csv", 
          row.names = FALSE)
write.csv(processed_data$log2, 
          "analysis_output/processed_log2_counts.csv", 
          row.names = FALSE)

###########################
# Print Analysis Summary
###########################

cat("\nAnalysis Summary:\n")
cat("Number of genes before filtering:", nrow(star_data), "\n")
cat("Number of genes after filtering:", nrow(processed_data$raw), "\n")
cat("Number of filtered genes:", processed_data$filtered_genes, "\n")
cat("Number of samples:", ncol(star_data) - 1, "\n")
cat("Number of significant genes (padj < 0.05):",
    sum(res_ordered$padj < 0.05, na.rm = TRUE), "\n")
cat("\nResults saved in 'analysis_output' directory\n")
