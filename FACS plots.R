# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Rsamtools", "GenomicAlignments", "GenomicFeatures"))

library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)

# Define file paths
annotation_file <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/Genes.gff3"
bam_dir <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/Unporechop demultiplexed fastq"
output_file <- "gene_counts.csv"

# List all BAM files in the directory
bam_files <- list.files(path = bam_dir, pattern = "*.bam$", full.names = TRUE)

# Create a TxDb object from GFF3 annotation file
txdb <- makeTxDbFromGFF(annotation_file, format = "gff3")
genes <- genes(txdb)

# Function to count reads
count_reads <- function(bam_file, genes) {
  param <- ScanBamParam(what = c("rname", "strand", "pos", "qwidth"))
  aln <- readGAlignments(bam_file, param = param)
  summarizeOverlaps(features = genes, reads = aln, mode = "Union", singleEnd = TRUE)
}

# Apply the function to all BAM files
count_list <- lapply(bam_files, count_reads, genes = genes)
count_matrix <- do.call(cbind, lapply(count_list, function(x) assays(x)$counts))
colnames(count_matrix) <- basename(bam_files)

# Save the count matrix to a CSV file
write.csv(as.data.frame(count_matrix), file = output_file, row.names = TRUE)

# Print completion message
cat("Gene counting completed. Results saved to:", output_file, "\n")
cat gene_counts.csv


# Install and load necessary packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ggplot2"))

library(ggplot2)

# Load gene counts from CSV file
gene_counts <- read.csv("gene_counts.csv", row.names = 1)

# Define groupings based on barcodes
group1 <- gene_counts[, 1:24]  # Barcodes 1-24 (uninfected-untreated)
group2 <- gene_counts[, 25:48] # Barcodes 25-48 (infected-untreated)
group3 <- gene_counts[, 49:72] # Barcodes 49-72 (infected-treated)
group4 <- gene_counts[, 73:92] # Barcodes 73-92 (uninfected-treated)

# Calculate mean expression for each group
group1_mean <- rowMeans(group1)
group2_mean <- rowMeans(group2)
group3_mean <- rowMeans(group3)
group4_mean <- rowMeans(group4)

# Prepare data frame for plotting
plot_data <- data.frame(
  GeneID = rownames(gene_counts),
  Group1 = group1_mean,
  Group2 = group2_mean,
  Group3 = group3_mean,
  Group4 = group4_mean
)

# Reshape data for plotting (tidy data format)
plot_data_long <- tidyr::gather(plot_data, key = "Group", value = "MeanExpression", -GeneID)

# Plot using ggplot2
ggplot(plot_data_long, aes(x = GeneID, y = MeanExpression, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean Gene Expression Across Groups",
       x = "Gene ID", y = "Mean Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Install and load necessary packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ggplot2", "pheatmap", "tidyr"))

library(ggplot2)
library(pheatmap)
library(tidyr)

# Load gene counts from CSV file
gene_counts <- read.csv("gene_counts.csv", row.names = 1)

# Define groupings based on barcodes
group1 <- gene_counts[, 1:24]  # Barcodes 1-24 (uninfected-untreated)
group2 <- gene_counts[, 25:48] # Barcodes 25-48 (infected-untreated)
group3 <- gene_counts[, 49:72] # Barcodes 49-72 (infected-treated)
group4 <- gene_counts[, 73:92] # Barcodes 73-92 (uninfected-treated)

# Calculate mean expression for each group
group1_mean <- rowMeans(group1)
group2_mean <- rowMeans(group2)
group3_mean <- rowMeans(group3)
group4_mean <- rowMeans(group4)

# Prepare data frame for plotting
plot_data <- data.frame(
  GeneID = rownames(gene_counts),
  Group1 = group1_mean,
  Group2 = group2_mean,
  Group3 = group3_mean,
  Group4 = group4_mean
)

# Reshape data for plotting (tidy data format)
plot_data_long <- tidyr::gather(plot_data, key = "Group", value = "MeanExpression", -GeneID)

# Bar Plot using ggplot2
ggplot(plot_data_long, aes(x = GeneID, y = MeanExpression, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean Gene Expression Across Groups",
       x = "Gene ID", y = "Mean Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Box Plot using ggplot2
# Reshape data to long format for individual samples
gene_counts_long <- tidyr::gather(gene_counts, key = "Barcode", value = "Expression")

# Add group labels based on barcodes
gene_counts_long$Group <- rep(c(rep("Uninfected-Untreated", 24),
                                rep("Infected-Untreated", 24),
                                rep("Infected-Treated", 24),
                                rep("Uninfected-Treated", 20)),
                              each = nrow(gene_counts))

# Boxplot
ggplot(gene_counts_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  labs(title = "Gene Expression Distribution Across Groups",
       x = "Group", y = "Expression")

# Heatmap using pheatmap
# Prepare data matrix for heatmap
data_matrix <- as.matrix(gene_counts)

# Generate heatmap
pheatmap(data_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         main = "Gene Expression Heatmap")


#Visualise FACS result
# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("flowCore")

install.packages("ggplot2")

library(flowCore)
library(ggplot2)

# Define the directory containing FCS files
fcs_directory <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/FACS"

# List all FCS files in the directory
fcs_files <- list.files(path = fcs_directory, pattern = "\\.fcs$", full.names = TRUE)

# Initialize an empty list to store data from all FCS files
all_fcs_data <- list()

# Loop through each FCS file
for (fcs_file in fcs_files) {
  # Read the FCS file
  fcs_data <- read.FCS(fcs_file)
  
  # Extract the expression data
  data <- exprs(fcs_data)
  
  # Add file identifier
  data <- as.data.frame(data)
  data$file_id <- basename(fcs_file)
  
  # Store the data in the list
  all_fcs_data[[basename(fcs_file)]] <- data
  
  # Print some basic info
  cat("Processed:", basename(fcs_file), "\n")
}

# Combine all data into a single data frame
combined_fcs_data <- do.call(rbind, all_fcs_data)

# Preview the combined data
head(combined_fcs_data)

# Plot FSC vs SSC for all files combined
ggplot(combined_fcs_data, aes(x = `FSC-A`, y = `SSC-A`, color = file_id)) +
  geom_point(alpha = 0.5) +
  labs(title = "FSC vs SSC Plot for All Files",
       x = "Forward Scatter (FSC-A)",
       y = "Side Scatter (SSC-A)") +
  theme_minimal()

# Create and save individual plots for each file
for (file_id in unique(combined_fcs_data$file_id)) {
  file_data <- subset(combined_fcs_data, file_id == file_id)
  
  plot <- ggplot(file_data, aes(x = `FSC-A`, y = `SSC-A`)) +
    geom_point(alpha = 0.5) +
    labs(title = paste("FSC vs SSC Plot for", file_id),
         x = "Forward Scatter (FSC-A)",
         y = "Side Scatter (SSC-A)") +
    theme_minimal()
  
  ggsave(filename = paste0("FSC_vs_SSC_", file_id, ".png"), plot = plot)
}
