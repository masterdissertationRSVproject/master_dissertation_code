library(ggplot2)

if (!file.exists("mapped_read_counts.txt")) {
  stop("Error: mapped_read_counts.txt file not found.")
}

read_counts <- read.table("mapped_read_counts.txt", header=FALSE, col.names=c("cell", "count"))

read_counts$cell_number <- as.numeric(gsub(".*_(\\d+)", "\\1", read_counts$cell))

# Debug
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

print("Data after condition assignment:")
print(head(read_counts))

print("Condition distribution:")
print(table(read_counts$condition))

read_counts <- read_counts[read_counts$condition != "Unknown", ]

print("Data after filtering unknown conditions:")
print(head(read_counts))

# Theme settings
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

# Remove existing file 
if (file.exists("dot_plot.png")) file.remove("dot_plot.png")

ggsave("dot_plot.png", plot=dot_plot)

# Box plot with log scale on y-axis
box_plot <- ggplot(read_counts, aes(x=condition, y=count, fill=condition)) +
  geom_boxplot() +
  scale_y_log10() +  # Log scale on y-axis
  common_theme +
  xlab("Condition") +
  ylab("Read Count") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file 
if (file.exists("box_plot.png")) file.remove("box_plot.png")

ggsave("box_plot.png", plot=box_plot)

# Violin plot with log scale on y-axis
violin_plot <- ggplot(read_counts, aes(x=condition, y=count, fill=condition)) +
  geom_violin() +
  scale_y_log10() +  # Log scale on y-axis
  common_theme +
  xlab("Condition") +
  ylab("Read Count") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file 
if (file.exists("violin_plot.png")) file.remove("violin_plot.png")

ggsave("violin_plot.png", plot=violin_plot)

# Histogram with log scale on y-axis
histogram <- ggplot(read_counts, aes(x=count, fill=condition)) +
  geom_histogram(bins=30, alpha=0.7, position="identity") +
  scale_x_log10() +  # Log scale on x-axis
  common_theme +
  xlab("Read Count") +
  ylab("Frequency") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file 
if (file.exists("histogram.png")) file.remove("histogram.png")

ggsave("histogram.png", plot=histogram)

# Density plot with log scale on x-axis
density_plot <- ggplot(read_counts, aes(x=count, fill=condition)) +
  geom_density(alpha=0.7) +
  scale_x_log10() +  # Log scale on x-axis
  common_theme +
  xlab("Read Count") +
  ylab("Density") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file 
if (file.exists("density_plot.png")) file.remove("density_plot.png")

ggsave("density_plot.png", plot=density_plot)

# Scatter plot
scatter_plot <- ggplot(read_counts, aes(x=cell_number, y=count, color=condition)) +
  geom_point() +
  scale_y_log10() +  # Log scale on y-axis
  common_theme +
  xlab("Cell Number") +
  ylab("Read Count") +
  ggtitle("RSV Read Counts per Condition")

# Remove existing file 
if (file.exists("scatter_plot.png")) file.remove("scatter_plot.png")

ggsave("scatter_plot.png", plot=scatter_plot)
