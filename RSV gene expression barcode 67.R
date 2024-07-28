# Load necessary libraries
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

# Set the paths to your files
bam_file <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/barcode_67.sorted_mapped.bam"
gff3_file <- "/Users/harintaka/nisrinafatiha/RSV Gene Expression/Genes.gff3"

# Create a transcript database from the GFF3 file
txdb <- makeTxDbFromGFF(gff3_file, format="gff3")

# Extract gene information from the database
genes <- genes(txdb)

# Read the BAM file
param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE))
bam_data <- readGAlignments(bam_file, param=param)

# Summarize overlaps to get read counts for each gene
se <- summarizeOverlaps(features=genes, reads=bam_data, mode="Union", ignore.strand=TRUE)

# Convert to a simple matrix of counts
count_matrix <- assay(se)

# Normalize the counts using CPM (Counts Per Million) for simplicity
cpm_matrix <- cpm(count_matrix)

# Create a dataframe with gene names and CPM values
gene_names <- rownames(count_matrix)
expression_data <- data.frame(Gene=gene_names, CPM=cpm_matrix)

# Print the expression data
print(expression_data)

# Visualize the gene expression profile
# Create a plot track for the genome
genome_track <- GenomeAxisTrack()

# Create a plot track for the gene regions
gene_track <- GeneRegionTrack(txdb, transcriptAnnotation="gene")

# Define the chromosome and coordinate range explicitly
chromosome <- "RSV_genome"  # Replace with the actual chromosome name
from <- 1
to <- max(end(genes))

# Create a plot track for the BAM alignments
alignments_track <- AlignmentsTrack(bam_file, isPaired=FALSE, chromosome=chromosome)

# Plot the tracks
plotTracks(list(genome_track, gene_track, alignments_track), from=from, to=to)

# Bar plot for selected genes
selected_genes <- c("F", "G", "L", "M", "N", "NS1", "NS2", "P", "SH") # Replace with actual gene names
selected_expression <- expression_data[expression_data$Gene %in% selected_genes, ]

ggplot(selected_expression, aes(x=Gene, y=reads)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  labs(title="Gene Expression Levels",
       x="Gene",
       y="Reads (CPM)")

# Heatmap for all genes
heatmap_data <- as.matrix(expression_data[, "reads"]) # Use only reads column
rownames(heatmap_data) <- expression_data$Gene

pheatmap(heatmap_data, 
         cluster_rows=TRUE, 
         cluster_cols=FALSE, 
         show_rownames=TRUE, 
         show_colnames=FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         main = "Gene Expression Heatmap")

#Trobleshoot 
# Install and load the rtracklayer package
if (!requireNamespace("rtracklayer", quietly = TRUE))
  BiocManager::install("rtracklayer")
library(rtracklayer)

# Import the GFF3 file
gff3_gr <- import.gff3(gff3_file)

# Inspect the imported GRanges object
head(gff3_gr)
# Inspect the GRanges object for potential issues
gff3_gr[seqnames(gff3_gr) == "gene-M2"]

# Filter out any problematic entries if necessary
gff3_filtered <- gff3_gr[!grepl("gene-M2", gff3_gr$ID)]

# Create TxDb from the filtered GRanges
txdb <- makeTxDbFromGRanges(gff3_filtered)

# Check the CDS entries that were dropped
cds_dropped <- subset(gff3_gr, ID %in% c("cds-ALS35591.1", "cds-ALS35592.1"))
cds_dropped

#Check the missing transcript
# Check all annotations related to gene-M2
gene_m2_annotations <- subset(gff3_gr, grepl("gene-M2", Parent) | grepl("gene-M2", ID))
gene_m2_annotations

# Create a new transcript entry for gene-M2
new_transcript <- GRanges(
  seqnames = Rle("KT992094.1"),
  ranges = IRanges(start = 7607, end = 8432), # Update the range as needed
  strand = Rle("+"),
  source = factor("Genbank"),
  type = factor("transcript"),
  score = NA_real_,
  phase = NA_integer_,
  ID = "transcript-M2",
  Name = "M2_transcript",
  gbkey = "mRNA",
  gene = "M2",
  gene_biotype = "protein_coding",
  Parent = "gene-M2",
  Dbxref = CharacterList(NA_character_),
  product = "M2 protein",
  protein_id = NA_character_
)

# Combine with the existing GRanges object
gff3_gr <- c(gff3_gr, new_transcript)

# Ensure the Parent field is correct for the CDS entries
mcols(gff3_gr)[which(mcols(gff3_gr)$ID %in% c("cds-ALS35591.1", "cds-ALS35592.1")), "Parent"] <- "transcript-M2"

# Create TxDb from the corrected GRanges object
txdb <- makeTxDbFromGRanges(gff3_gr)

# Inspect the current GRanges object
gff3_gr

# Check for exons related to gene-M2
gene_m2_exons <- subset(gff3_gr, grepl("transcript-M2", Parent) | grepl("transcript-M2", ID))
gene_m2_exons

# If exons are missing, manually add them
new_exons <- GRanges(
  seqnames = Rle("KT992094.1"),
  ranges = IRanges(start = c(7607, 8160), end = c(8191, 8432)), # Update the ranges as needed
  strand = Rle("+"),
  source = factor("Genbank"),
  type = factor("exon"),
  score = NA_real_,
  phase = NA_integer_,
  ID = c("exon-M2-1", "exon-M2-2"),
  Name = c("exon1", "exon2"),
  gbkey = "exon",
  gene = "M2",
  gene_biotype = "protein_coding",
  Parent = "transcript-M2",
  Dbxref = CharacterList(NA_character_),
  product = NA_character_,
  protein_id = NA_character_
)

# Combine with the existing GRanges object
gff3_gr <- c(gff3_gr, new_exons)

# Ensure the Parent field is correct for the CDS entries
mcols(gff3_gr)[which(mcols(gff3_gr)$ID %in% c("cds-ALS35591.1", "cds-ALS35592.1")), "Parent"] <- "transcript-M2"

# Create TxDb from the corrected GRanges object
txdb <- makeTxDbFromGRanges(gff3_gr)

# Inspect the current GRanges object for all annotations related to transcript-M2
gene_m2_annotations <- subset(gff3_gr, grepl("transcript-M2", Parent) | grepl("transcript-M2", ID))
gene_m2_annotations

# Create exon entries for transcript-M2 if not already present
new_exons <- GRanges(
  seqnames = Rle("KT992094.1"),
  ranges = IRanges(start = c(7607, 8160), end = c(8191, 8432)), # Update the ranges as needed
  strand = Rle("+"),
  source = factor("Genbank"),
  type = factor("exon"),
  score = NA_real_,
  phase = NA_integer_,
  ID = c("exon-M2-1", "exon-M2-2"),
  Name = c("exon1", "exon2"),
  gbkey = "exon",
  gene = "M2",
  gene_biotype = "protein_coding",
  Parent = "transcript-M2",
  Dbxref = CharacterList(NA_character_),
  product = NA_character_,
  protein_id = NA_character_
)

# Combine with the existing GRanges object
gff3_gr <- c(gff3_gr, new_exons)

# Ensure the Parent field is correct for the CDS entries
mcols(gff3_gr)[which(mcols(gff3_gr)$ID %in% c("cds-ALS35591.1", "cds-ALS35592.1")), "Parent"] <- "transcript-M2"

# Check the updated GRanges object
gff3_gr

# Create TxDb from the corrected GRanges object
txdb <- makeTxDbFromGRanges(gff3_gr)

# Verify that exons are linked correctly
subset(gff3_gr, type == "exon" & Parent == "transcript-M2")

# Verify that CDS entries are linked correctly
subset(gff3_gr, type == "CDS" & Parent == "transcript-M2")

# Verify the transcript entry
subset(gff3_gr, type == "transcript" & ID == "transcript-M2")

