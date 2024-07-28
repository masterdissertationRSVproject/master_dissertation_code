# Load necessary library
library(ggplot2)

# Read data from the file
read_counts <- read.table("mapped_read_counts.txt", header=FALSE, col.names=c("barcode", "count"))

# Extract numeric part of barcode identifiers
read_counts$cell_number <- as.numeric(gsub("barcode_(\\d+)\\.sorted", "\\1", read_counts$barcode))

# Debug: Identify any NAs introduced during coercion
na_cells <- read_counts[is.na(read_counts$cell_number), "barcode"]
if (length(na_cells) > 0) {
  print("Warning: NAs introduced by coercion. Problematic cells:")
  print(na_cells)
}

# Assign conditions based on cell numbers
read_counts$condition <- ifelse(read_counts$cell_number %in% 1:24, "Uninfected Untreated",
                                ifelse(read_counts$cell_number %in% 25:48, "Infected Untreated",
                                       ifelse(read_counts$cell_number %in% 49:72, "Infected Treated",
                                              ifelse(read_counts$cell_number %in% 73:96, "Uninfected Treated",
                                                     "Unknown"))))

# Filter out "Unknown" condition (optional, if needed)
read_counts <- read_counts[read_counts$condition != "Unknown", ]

# Create violin plot
violin_plot <- ggplot(read_counts, aes(x=condition, y=count, fill=condition)) +
  geom_violin() +
  scale_y_log10() +  # Log scale on y-axis
  theme_minimal() +
  xlab("Condition") +
  ylab("Read Count (log10)") +
  ggtitle("Violin Plot of Read Counts by Condition")

# Save the plot
ggsave("violin_plot.png", plot=violin_plot)

# Display the plot
print(violin_plot)
