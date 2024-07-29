# Load libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "GenomicAlignments", "Gviz", "Rsamtools", "edgeR", "pheatmap", "RColorBrewer", "ggplot2"))

library(GenomicFeatures)
library(GenomicAlignments)
library(Gviz)
library(Rsamtools)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# Set the paths 
bam_file <- "path_to_my_file/barcode_67.sorted_mapped.bam"
gff3_file <- "path_to_annotation_file/Genes.gff3"

# Create a transcript database from the GFF3 file
txdb <- makeTxDbFromGFF(gff3_file, format="gff3")

# Extract gene information from the database
genes <- genes(txdb)

# Read the BAM file
param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE))
bam_data <- readGAlignments(bam_file, param=param)

# Summarise overlaps to get read counts for each gene
se <- summarizeOverlaps(features=genes, reads=bam_data, mode="Union", ignore.strand=TRUE)

# Convert to a simple matrix of counts
count_matrix <- assay(se)

# Normalise the counts using CPM (Counts Per Million) 
cpm_matrix <- cpm(count_matrix)

# Create a dataframe with gene names and CPM values
gene_names <- rownames(count_matrix)
expression_data <- data.frame(Gene=gene_names, CPM=cpm_matrix)

# Print the expression data
print(expression_data)

# Visualise the gene expression profile
# Create a plot track for the genome
genome_track <- GenomeAxisTrack()

# Create a plot track for the gene regions
gene_track <- GeneRegionTrack(txdb, transcriptAnnotation="gene")

# Define the chromosome and coordinate range explicitly
chromosome <- "RSV_genome"  
from <- 1
to <- max(end(genes))

# Create a plot track for the BAM alignments
alignments_track <- AlignmentsTrack(bam_file, isPaired=FALSE, chromosome=chromosome)

# Plot tracks
plotTracks(list(genome_track, gene_track, alignments_track), from=from, to=to)

# Bar plot for selected genes
selected_genes <- c("F", "G", "L", "M", "N", "NS1", "NS2", "P", "SH") 
selected_expression <- expression_data[expression_data$Gene %in% selected_genes, ]

ggplot(selected_expression, aes(x=Gene, y=reads)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  labs(title="Gene Expression Levels",
       x="Gene",
       y="Reads (CPM)")

# Heatmap for all genes
heatmap_data <- as.matrix(expression_data[, "reads"]) 
rownames(heatmap_data) <- expression_data$Gene

pheatmap(heatmap_data, 
         cluster_rows=TRUE, 
         cluster_cols=FALSE, 
         show_rownames=TRUE, 
         show_colnames=FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         main = "Gene Expression Heatmap")
