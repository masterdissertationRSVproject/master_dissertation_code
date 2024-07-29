library(ggplot2)

read_counts <- read.table("mapped_read_counts.txt", header=FALSE, col.names=c("barcode", "count"))

read_counts$cell_number <- as.numeric(gsub(".*_(\\d+)", "\\1", read_counts$barcode))

# Assign conditions based on cell numbers
read_counts$condition <- ifelse(read_counts$cell_number %in% 1:24, "Uninfected",
                                ifelse(read_counts$cell_number %in% 25:48, "RSV",
                                       ifelse(read_counts$cell_number %in% 49:72, "RSV+DAP",
                                              ifelse(read_counts$cell_number %in% 73:92, "DAP",
                                                     "Unknown"))))

read_counts <- read_counts[read_counts$condition != "Unknown", ]

# Create violin plot
violin_plot <- ggplot(read_counts, aes(x=condition, y=count, fill=condition)) +
  geom_violin() +
  scale_y_log10() +  # Log scale on y-axis
  theme_minimal() +
  xlab("Condition") +
  ylab("Read Count (log10)") +
  ggtitle("RSV Read Counts by Condition")

# Save the plot
ggsave("violin_plot.png", plot=violin_plot)

# Display the plot
print(violin_plot)
