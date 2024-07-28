# Load necessary library
library(ggplot2)

# Set up graphics device
png()  # Or pdf(), depending on your preferred format

# Check if mapped_read_counts.txt file exists
if (!file.exists("mapped_read_counts.txt")) {
  stop("Error: mapped_read_counts.txt file not found.")
}

# Load read counts from file
read_counts <- read.table("mapped_read_counts.txt", header=FALSE, col.names=c("cell", "count"))

# Debug: Print the first few lines of read_counts to verify data loading
print("Initial data loaded:")
print(head(read_counts))

# Extract numeric part of cell identifiers
read_counts$cell_number <- as.numeric(gsub(".*_(\\d+)", "\\1", read_counts$cell))

# Debug: Identify any NAs introduced during coercion
na_cells <- read_counts[is.na(read_counts$cell_number), "cell"]
if (length(na_cells) > 0) {
  print("Warning: NAs introduced by coercion. Problematic cells:")
  print(na_cells)
}

# Assign conditions based on cell numbers
read_counts$condition <- ifelse(read_counts$cell_number %in% 1:24, "Uninfected",
                                ifelse(read_counts$cell_number %in% 25:48, "RSV",
                                       ifelse(read_counts$cell_number %in% 49:72, "RSV+DAP",
                                              ifelse(read_counts$cell_number %in% 73:96, "DAP",
                                                     "Unknown"))))

# Debug: Print the first few lines of read_counts to verify condition assignment
print("Data after condition assignment:")
print(head(read_counts))

# Debug: Check the distribution of conditions
print("Condition distribution:")
print(table(read_counts$condition))

# Filter out "Unknown" condition (optional, if needed)
read_counts <- read_counts[read_counts$condition != "Unknown", ]

# Debug: Print the first few lines of read_counts to verify filtering
print("Data after filtering unknown conditions:")
print(head(read_counts))

# Common theme settings
common_theme <- theme_minimal() +
  theme(text = element_text(family = "Arial", size = 12))

# Dot plot with log scale on y-axis
dot_plot <- ggplot(read_counts, aes(x=condition, y=count, color=condition)) +
  geom_point() +
  scale_y_log10() +  # Log scale on y-axis
  common_theme +
  xlab("Condition") +
  ylab("Read Count") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file if present
if (file.exists("dot_plot.png")) file.remove("dot_plot.png")

# Save dot plot
ggsave("dot_plot.png", plot=dot_plot)

# Box plot with log scale on y-axis
box_plot <- ggplot(read_counts, aes(x=condition, y=count, fill=condition)) +
  geom_boxplot() +
  scale_y_log10() +  # Log scale on y-axis
  common_theme +
  xlab("Condition") +
  ylab("Read Count") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file if present
if (file.exists("box_plot.png")) file.remove("box_plot.png")

# Save box plot
ggsave("box_plot.png", plot=box_plot)

# Violin plot with log scale on y-axis
violin_plot <- ggplot(read_counts, aes(x=condition, y=count, fill=condition)) +
  geom_violin() +
  scale_y_log10() +  # Log scale on y-axis
  common_theme +
  xlab("Condition") +
  ylab("Read Count") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file if present
if (file.exists("violin_plot.png")) file.remove("violin_plot.png")

# Save violin plot
ggsave("violin_plot.png", plot=violin_plot)

# Histogram with log scale on y-axis
histogram <- ggplot(read_counts, aes(x=count, fill=condition)) +
  geom_histogram(bins=30, alpha=0.7, position="identity") +
  scale_x_log10() +  # Log scale on x-axis
  common_theme +
  xlab("Read Count") +
  ylab("Frequency") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file if present
if (file.exists("histogram.png")) file.remove("histogram.png")

# Save histogram
ggsave("histogram.png", plot=histogram)

# Density plot with log scale on x-axis
density_plot <- ggplot(read_counts, aes(x=count, fill=condition)) +
  geom_density(alpha=0.7) +
  scale_x_log10() +  # Log scale on x-axis
  common_theme +
  xlab("Read Count") +
  ylab("Density") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file if present
if (file.exists("density_plot.png")) file.remove("density_plot.png")

# Save density plot
ggsave("density_plot.png", plot=density_plot)

# Scatter plot
scatter_plot <- ggplot(read_counts, aes(x=cell_number, y=count, color=condition)) +
  geom_point() +
  scale_y_log10() +  # Log scale on y-axis
  common_theme +
  xlab("Cell Number") +
  ylab("Read Count") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file if present
if (file.exists("scatter_plot.png")) file.remove("scatter_plot.png")

# Save scatter plot
ggsave("scatter_plot.png", plot=scatter_plot)

# Close graphics device
dev.off()

print("All plots generated successfully.")
