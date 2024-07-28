# Load necessary library
library(ggplot2)

# Create the data frame
data <- data.frame(
  gene_name = c("NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L"),
  gene_length = c(532, 503, 1203, 914, 958, 410, 923, 1903, 961, 6578),
  coverage = c(303, 538, 1109, 978, 615, 327, 1278, 640, 426, 52)
)

# Print the data frame to ensure it is created correctly
print(data)

# Create the dot plot with different colors for each gene
plot <- ggplot(data, aes(x = coverage, y = gene_length, color = gene_name)) +
  geom_point(size = 3) +  # Set size of points
  labs(
    x = "Coverage (total count)",
    y = "Gene Length (bp)",
    title = "Dot Plot of RSV Gene Length vs. Coverage"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("NS1" = "red", "NS2" = "blue", "N" = "green",
                                "P" = "purple", "M" = "orange", "SH" = "pink",
                                "G" = "cyan", "F" = "magenta", "M2" = "brown",
                                "L" = "black"))

# Print the plot
print(plot)


# Load necessary library
library(ggplot2)

# Create the data frame
data <- data.frame(
  gene_name = c("NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L"),
  gene_length = c(532, 503, 1203, 914, 958, 410, 923, 1903, 961, 6578),
  coverage = c(303, 538, 1109, 978, 615, 327, 1278, 640, 426, 52)
)

# Print the data frame to ensure it is created correctly
print(data)

# Create the scatter plot with a regression line and correlation coefficient
plot <- ggplot(data, aes(x = coverage, y = gene_length, color = gene_name)) +
  geom_point(size = 3) +  # Set size of points
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Add regression line
  labs(
    x = "Coverage (total count)",
    y = "Gene Length (bp)",
    title = "Scatter Plot of RSV Gene Length vs Coverage"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("NS1" = "red", "NS2" = "blue", "N" = "green",
                                "P" = "purple", "M" = "orange", "SH" = "pink",
                                "G" = "cyan", "F" = "magenta", "M2" = "brown",
                                "L" = "black"))

# Print the plot
print(plot)

# Calculate and print correlation coefficient
correlation <- cor(data$coverage, data$gene_length)
print(paste("Correlation coefficient:", round(correlation, 2)))


# Load necessary library
library(ggplot2)

# Create the data frame
data <- data.frame(
  gene_name = c("NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L"),
  gene_length = c(532, 503, 1203, 914, 958, 410, 923, 1903, 961, 6578),
  coverage = c(303, 538, 1109, 978, 615, 327, 1278, 640, 426, 52)
)

# Print the data frame to ensure it is created correctly
print(data)

# Create the scatter plot with a regression line and correlation coefficient
plot <- ggplot(data, aes(x = gene_length, y = coverage, color = gene_name)) +
  geom_point(size = 3) +  # Set size of points
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Add regression line
  labs(
    x = "Gene Length (bp)",
    y = "Coverage (total count)",
    title = "Scatter Plot of RSV Gene Length vs Coverage"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("NS1" = "red", "NS2" = "blue", "N" = "green",
                                "P" = "purple", "M" = "orange", "SH" = "pink",
                                "G" = "cyan", "F" = "magenta", "M2" = "brown",
                                "L" = "black"))

# Print the plot
print(plot)

# Calculate and print correlation coefficient
correlation <- cor(data$gene_length, data$coverage)
print(paste("Correlation coefficient:", round(correlation, 2)))
