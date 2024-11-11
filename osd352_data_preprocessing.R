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

###############################
# Data Loading & Pre-processing
###############################

# Function to load and clean STAR output
preprocess_rnaseq <- function(star_data) {
  gene_id_col <- colnames(star_data)[1]
  expressed_genes <- star_data[rowSums(star_data[,-1] > 0) >= (ncol(star_data)-1)*0.25, ]
  log2_counts <- log2(expressed_genes[,-1] + 1)
  log2_counts[[gene_id_col]] <- expressed_genes[[gene_id_col]]
  return(list(
    raw = expressed_genes,
    log2 = log2_counts,
    gene_id_col = gene_id_col
  ))
}

# Function to plot expression density
plot_expression_density <- function(log2_counts, gene_id_col) {
  melted_counts <- log2_counts %>%
    select(-all_of(gene_id_col)) %>%
    gather(key = "Sample", value = "Log2_Counts")
  
  ggplot(melted_counts, aes(x = Log2_Counts, color = Sample)) +
    geom_density() +
    theme_minimal() +
    labs(title = "Distribution of Log2 Expression Values",
         x = "Log2(Counts + 1)",
         y = "Density") +
    theme(legend.position = "none")
}

# Function to plot sample correlation heatmap
plot_sample_correlation <- function(log2_counts, gene_id_col) {
  cor_matrix <- cor(log2_counts %>% select(-all_of(gene_id_col)))
  pheatmap(cor_matrix,
           main = "Sample Correlation Heatmap",
           show_rownames = TRUE,
           show_colnames = TRUE,
           display_numbers = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(50))
}

# Function to plot gene statistics
plot_gene_stats <- function(processed_data) {
  gene_stats <- data.frame(
    Mean_Expression = rowMeans(processed_data$log2 %>% select(-all_of(processed_data$gene_id_col))),
    CV = apply(processed_data$log2 %>% select(-all_of(processed_data$gene_id_col)), 1, sd) / 
      apply(processed_data$log2 %>% select(-all_of(processed_data$gene_id_col)), 1, mean)
  )
  ggplot(gene_stats, aes(x = Mean_Expression, y = CV)) +
    geom_point(alpha = 0.2) +
    theme_minimal() +
    labs(title = "Gene Expression Variability",
         x = "Mean Log2 Expression",
         y = "Coefficient of Variation")
}

# Setting dynamic working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading the data
star_data <- read.csv("./data/GLDS-352_rna_seq_STAR_Unnormalized_Counts.csv")

# Processing the data
processed_data <- preprocess_rnaseq(star_data)

# Creating the plots
p1 <- plot_expression_density(processed_data$log2, processed_data$gene_id_col)
p2 <- plot_gene_stats(processed_data)
plot_sample_correlation(processed_data$log2, processed_data$gene_id_col)

# Some additional information
cat("Number of genes before filtering:", nrow(star_data), "\n")
cat("Number of genes after filtering:", nrow(processed_data$raw), "\n")
cat("Number of samples:", ncol(star_data)-1, "\n")

###########################
# Saving the processed data
###########################

write.csv(processed_data$raw, "./data/processed_counts.csv", row.names = FALSE)
write.csv(processed_data$log2, "./data/processed_log2_counts.csv", row.names = FALSE)
