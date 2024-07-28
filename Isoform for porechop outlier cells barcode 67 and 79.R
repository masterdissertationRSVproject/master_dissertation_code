# Load necessary libraries
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(data.table)
library(ggplot2)

# Define file paths (update these paths according to your actual files)
bam_file <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/Unporechop demultiplexed fastq/barcode_67.sorted_mapped.bam"
bed_file <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/Unporechop demultiplexed fastq/barcode_67.sorted_mapped.bed"
genes_bed_file <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/Genes.gff3"
intersected_reads_file <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/intersected_reads.bed"
counts_per_isoform_file <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/counts_per_isoform.bed"

# Ensure Bedtools is available
bedtools_path <- system("which bedtools", intern = TRUE)
if (bedtools_path == "") {
  stop("Bedtools is not installed or not in the PATH. Please install Bedtools and try again.")
}

# Convert BAM to BED
cmd1 <- paste0(bedtools_path, " bamtobed -i '", bam_file, "' > '", bed_file, "'")
system(cmd1)

# Intersect reads with annotations
cmd2 <- paste0(bedtools_path, " intersect -a '", bed_file, "' -b '", genes_bed_file, "' > '", intersected_reads_file, "'")
system(cmd2)

# Summarize intersections
cmd3 <- paste0(bedtools_path, " intersect -a '", genes_bed_file, "' -b '", bed_file, "' -c > '", counts_per_isoform_file, "'")
system(cmd3)

# Load and inspect the intersected reads
intersected_reads <- fread(intersected_reads_file)
print(head(intersected_reads))

# Load and inspect the counts per isoform
counts_per_isoform <- fread(counts_per_isoform_file)
print(head(counts_per_isoform))

# Ensure the V10 column has appropriate numeric values
counts_per_isoform$V10 <- as.numeric(counts_per_isoform$V10)
summary(counts_per_isoform$V10)

# Plot the distribution of read counts per gene as a dot plot
ggplot(counts_per_isoform, aes(x = V10)) +
  geom_dotplot(binwidth = 1, fill = "blue", color = "black", dotsize = 0.5) +
  labs(title = "Distribution of Read Counts per Gene",
       x = "Read Count",
       y = "Frequency") +
  theme_minimal()
