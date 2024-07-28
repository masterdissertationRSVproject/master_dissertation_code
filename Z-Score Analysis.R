# Load necessary libraries
library(ggplot2)
library(outliers)
library(MASS)

setwd("/Users/harintaka/nisrinafatiha/RSV Gene Expression")

# Your read counts data
read_counts <- data.frame(
  cell_number = 1:96,
  count = c(234, 35, 52, 17, 42, 40, 70, 20, 50, 54, 46, 52, 15, 17, 32, 46, 39, 63, 13, 25, 18, 44, 24, 41, 34, 48, 34, 25, 30, 16, 6, 61, 59, 12, 130, 17, 36, 27, 54, 41, 69, 19, 16, 51, 84, 14, 12, 24, 123, 14, 56, 17, 35, 45, 20, 42, 17, 23, 18, 39, 0, 88, 62, 148, 153, 146, 77356, 81, 5, 29, 25, 18, 51, 71, 49, 27, 1027, 210, 51595, 21, 18, 16, 7, 69, 44, 155, 117, 88, 144, 146, 5, 38, 0, 0, 0, 0)
)

# Plotting the read counts
ggplot(read_counts, aes(x = cell_number, y = count)) +
  geom_point() +
  scale_y_log10() +
  theme_minimal() +
  xlab("Cell Number") +
  ylab("Read Count (log10)") +
  ggtitle("Read Counts per Cell")

# Z-Score Calculation
read_counts$z_score <- (read_counts$count - mean(read_counts$count)) / sd(read_counts$count)

# Outlier Detection: Grubbs' Test
grubbs_result <- grubbs.test(read_counts$count)

# Outlier Detection: Modified Z-Score
threshold <- 3.5
read_counts$modified_z <- 0.6745 * (read_counts$count - median(read_counts$count)) / mad(read_counts$count)
outliers <- read_counts[abs(read_counts$modified_z) > threshold, ]

# Hypothesis Testing: Mann-Whitney U Test
# Separate high and low reads for testing
high_reads <- read_counts$count[read_counts$count > 1000]
low_reads <- read_counts$count[read_counts$count <= 1000]
test_result <- wilcox.test(high_reads, low_reads)

# Background Noise Modeling: Negative Binomial
# Fit a negative binomial distribution to the low reads data
fit <- fitdistr(low_reads, "negative binomial")

# Print results
print("Grubbs' Test Result:")
print(grubbs_result)

print("Outliers Based on Modified Z-Score:")
print(outliers)

print("Mann-Whitney U Test Result:")
print(test_result)

print("Negative Binomial Fit Summary:")
print(summary(fit))


#Insert the uninfected and infected into account 
# Load necessary libraries
library(MASS)
library(ggplot2)

# Read counts data (replace this with your actual data loading method)
read_counts <- data.frame(
  cell_number = 1:96,
  count = c(
    0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 2, 1, 1, 3, 1, 1, 0, 5, 1, 0, 0, 1, 0, 2, 2, 0, 0, 4, 0, 1, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 1, 0, 10, 1, 0, 5872, 0, 1, 0, 2, 0, 1, 1, 0, 0, 75, 16, 4099, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 5, 0, 3, 0, 0, 0, 0
  )
)

# Create a vector to indicate whether a barcode is uninfected or not
read_counts$infected <- ifelse(read_counts$cell_number %in% c(1:24, 73:96), "Uninfected", "Infected")

# Separate the data
uninfected_counts <- read_counts[read_counts$infected == "Uninfected", ]
infected_counts <- read_counts[read_counts$infected == "Infected", ]

# 1. Quality Control and Background Noise Modeling
# Fit a Negative Binomial model to uninfected counts
fit_background <- fitdistr(uninfected_counts$count, "negative binomial")

# Print fit results
print("Negative Binomial Fit Summary:")
print(summary(fit_background))

# 2. Hypothesis Testing
# Extract the fitted distribution parameters
size <- fit_background$estimate["size"]
mu <- fit_background$estimate["mu"]

# Generate the expected background distribution for comparison
expected_background <- dnbinom(infected_counts$count, size = size, mu = mu)

# Perform Kolmogorov-Smirnov Test
ks_test <- ks.test(infected_counts$count, expected_background)
print("Kolmogorov-Smirnov Test Result:")
print(ks_test)

# 3. Visual Inspection
# Plot histogram of read counts for infected and uninfected cells
ggplot(read_counts, aes(x = count, fill = infected)) +
  geom_histogram(binwidth = 50, position = "dodge") +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Read Counts Distribution",
       x = "Read Count",
       y = "Frequency")

# 4. Statistical Analysis of Outliers
# Calculate Z-scores for infected counts
infected_counts$z_score <- (infected_counts$count - mean(uninfected_counts$count)) / sd(uninfected_counts$count)

# Identify extreme Z-scores
extreme_z <- infected_counts[infected_counts$z_score > 3, ]
print("Extreme Z-Scores:")
print(extreme_z)

# Calculate Modified Z-scores
infected_counts$modified_z <- 0.6745 * (infected_counts$count - median(uninfected_counts$count)) / mad(uninfected_counts$count)
threshold <- 3.5
outliers <- infected_counts[abs(infected_counts$modified_z) > threshold, ]
print("Outliers Based on Modified Z-Scores:")
print(outliers)

# Mann-Whitney U Test
high_reads <- infected_counts$count[infected_counts$count > 1000]
low_reads <- infected_counts$count[infected_counts$count <= 1000]
test_result <- wilcox.test(high_reads, low_reads)
print("Mann-Whitney U Test Result:")
print(test_result)

# Print results of the analysis
print("Results Summary:")
print(list(
  "Negative Binomial Fit" = summary(fit_background),
  "Kolmogorov-Smirnov Test" = ks_test,
  "Extreme Z-Scores" = extreme_z,
  "Outliers" = outliers,
  "Mann-Whitney U Test" = test_result
))

#codes incorporating the 0 value!!!!!
# Load necessary libraries
library(MASS)
library(ggplot2)
library(pscl)  # For zero-inflated models if needed

# Sample data
read_counts <- data.frame(
  cell_number = 1:96,
  count = c(
0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 2, 1, 1, 3, 1, 1, 0, 5, 1, 0, 0, 1, 0, 2, 2, 0, 0, 4, 0, 1, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 1, 0, 10, 1, 0, 5872, 0, 1, 0, 2, 0, 1, 1, 0, 0, 75, 16, 4099, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 5, 0, 3, 0, 0, 0, 0
)
)
# Create a vector to indicate whether a barcode is uninfected or not
read_counts$infected <- ifelse(read_counts$cell_number %in% c(1:24, 73:96), "Uninfected", "Infected")

# Separate the data
uninfected_counts <- read_counts[read_counts$infected == "Uninfected", ]
infected_counts <- read_counts[read_counts$infected == "Infected", ]

# Add a small constant to avoid zero counts
uninfected_counts$count <- uninfected_counts$count + 1

# Fit a Negative Binomial model to filtered uninfected counts
fit_background <- tryCatch(
  fitdistr(uninfected_counts$count, "negative binomial"),
  error = function(e) NULL
)

# If fitting fails, use zero-inflated negative binomial as an alternative
if (is.null(fit_background)) {
  fit_background <- tryCatch(
    zeroinfl(count ~ 1 | 1, data = uninfected_counts, dist = "negbin"),
    error = function(e) NULL
  )
}

# Background Noise Analysis
if (!is.null(fit_background)) {
  print(summary(fit_background))
}

# Plot histogram of counts
hist(uninfected_counts$count, breaks = 30, main = "Histogram of Uninfected Counts", xlab = "Count")

# Plot histogram of infected counts
hist(infected_counts$count, breaks = 30, main = "Histogram of Infected Counts", xlab = "Count")

# Kolmogorov-Smirnov Test
if (!is.null(fit_background) && "negbin" %in% class(fit_background)) {
  size <- fit_background$estimate["size"]
  mu <- fit_background$estimate["mu"]
  expected_background <- dnbinom(infected_counts$count, size = size, mu = mu)
  
  ks_test <- ks.test(infected_counts$count, expected_background)
  print(ks_test)
}

# Mann-Whitney U Test
high_reads <- infected_counts$count[infected_counts$count > 1000]
low_reads <- infected_counts$count[infected_counts$count <= 1000]

test_result <- wilcox.test(high_reads, low_reads)
print(test_result)

# Plot distribution of read counts
ggplot(read_counts, aes(x = count, fill = infected)) +
  geom_histogram(binwidth = 50, position = "dodge") +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Read Counts Distribution",
       x = "Read Count",
       y = "Frequency")

# Z-Score Analysis
mean_uninfected <- mean(uninfected_counts$count)
sd_uninfected <- sd(uninfected_counts$count)
infected_counts$z_score <- (infected_counts$count - mean_uninfected) / sd_uninfected

# Identify extreme Z-scores
extreme_z <- infected_counts[infected_counts$z_score > 3, ]
print(extreme_z)

#Visualisation of the data above 
# Load necessary libraries
install.packages(c("ggplot2", "MASS", "pscl", "outliers"))
library(ggplot2)
library(MASS)
library(pscl)
library(outliers)

# Prepare the data
read_counts <- data.frame(
  cell_number = 1:96,
  count = c(234, 35, 52, 17, 42, 40, 70, 20, 50, 54, 46, 52, 15, 17, 32, 46, 39, 63, 13, 25, 18, 44, 24, 41, 34, 48, 34, 25, 30, 16, 6, 61, 59, 12, 130, 17, 36, 27, 54, 41, 69, 19, 16, 51, 84, 14, 12, 24, 123, 14, 56, 17, 35, 45, 20, 42, 17, 23, 18, 39, 0, 88, 62, 148, 153, 146, 77356, 81, 5, 29, 25, 18, 51, 71, 49, 27, 1027, 210, 51595, 21, 18, 16, 7, 69, 44, 155, 117, 88, 144, 146, 5, 38, 0, 0, 0, 0)
)

read_counts$infected <- ifelse(read_counts$cell_number %in% c(1:24, 73:96), "Uninfected", "Infected")

uninfected_counts <- read_counts[read_counts$infected == "Uninfected", ]
infected_counts <- read_counts[read_counts$infected == "Infected", ]

# Fit a Zero-Inflated Negative Binomial model to uninfected counts
fit_background <- tryCatch(
  zeroinfl(count ~ 1 | 1, data = uninfected_counts, dist = "negbin"),
  error = function(e) NULL
)

# Plot histograms
ggplot(read_counts, aes(x = count, fill = infected)) +
  geom_histogram(binwidth = 50, position = "dodge", alpha = 0.7) +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Read Counts Distribution",
       x = "Read Count",
       y = "Frequency",
       fill = "Infection Status")

# Kolmogorov-Smirnov Test
if (!is.null(fit_background) && "negbin" %in% class(fit_background)) {
  size <- fit_background$theta
  mu <- fit_background$coefficients$count
  
  expected_background <- dnbinom(infected_counts$count, size = size, mu = mu)
  ks_test <- ks.test(infected_counts$count, expected_background)
  print(ks_test)
}

# Mann-Whitney U Test
high_reads <- infected_counts$count[infected_counts$count > 1000]
low_reads <- infected_counts$count[infected_counts$count <= 1000]
test_result <- wilcox.test(high_reads, low_reads)
print(test_result)

# Z-Score Analysis
mean_uninfected <- mean(uninfected_counts$count)
sd_uninfected <- sd(uninfected_counts$count)

infected_counts$z_score <- (infected_counts$count - mean_uninfected) / sd_uninfected
extreme_z <- infected_counts[infected_counts$z_score > 3, ]
print(extreme_z)

# Additional Plots for Dissertation
# Load necessary libraries
install.packages(c("ggplot2", "MASS", "pscl", "outliers"))
library(ggplot2)
library(MASS)
library(pscl)
library(outliers)

# Prepare the data
read_counts <- data.frame(
  cell_number = 1:96,
  count = c(234, 35, 52, 17, 42, 40, 70, 20, 50, 54, 46, 52, 15, 17, 32, 46, 39, 63, 13, 25, 18, 44, 24, 41, 34, 48, 34, 25, 30, 16, 6, 61, 59, 12, 130, 17, 36, 27, 54, 41, 69, 19, 16, 51, 84, 14, 12, 24, 123, 14, 56, 17, 35, 45, 20, 42, 17, 23, 18, 39, 0, 88, 62, 148, 153, 146, 77356, 81, 5, 29, 25, 18, 51, 71, 49, 27, 1027, 210, 51595, 21, 18, 16, 7, 69, 44, 155, 117, 88, 144, 146, 5, 38, 0, 0, 0, 0)
)

read_counts$infected <- ifelse(read_counts$cell_number %in% c(1:24, 73:96), "Uninfected", "Infected")

uninfected_counts <- read_counts[read_counts$infected == "Uninfected", ]
infected_counts <- read_counts[read_counts$infected == "Infected", ]

# Fit a Zero-Inflated Negative Binomial model to uninfected counts
fit_background <- tryCatch(
  zeroinfl(count ~ 1 | 1, data = uninfected_counts, dist = "negbin"),
  error = function(e) NULL
)

# Plot histograms
plot1 <- ggplot(read_counts, aes(x = count, fill = infected)) +
  geom_histogram(binwidth = 50, position = "dodge", alpha = 0.7) +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Read Counts Distribution",
       x = "Read Count",
       y = "Frequency",
       fill = "Infection Status")
print(plot1)

# Kolmogorov-Smirnov Test
if (!is.null(fit_background) && "negbin" %in% class(fit_background)) {
  size <- fit_background$theta
  mu <- fit_background$coefficients$count
  
  expected_background <- dnbinom(infected_counts$count, size = size, mu = mu)
  ks_test <- ks.test(infected_counts$count, expected_background)
  print(ks_test)
}

# Mann-Whitney U Test
high_reads <- infected_counts$count[infected_counts$count > 1000]
low_reads <- infected_counts$count[infected_counts$count <= 1000]
test_result <- wilcox.test(high_reads, low_reads)
print(test_result)

# Z-Score Analysis
mean_uninfected <- mean(uninfected_counts$count)
sd_uninfected <- sd(uninfected_counts$count)

infected_counts$z_score <- (infected_counts$count - mean_uninfected) / sd_uninfected
extreme_z <- infected_counts[infected_counts$z_score > 3, ]
print(extreme_z)

# Additional Plots for Dissertation
plot2 <- ggplot(read_counts, aes(x = infected, y = count, fill = infected)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Boxplot of Read Counts by Infection Status",
       x = "Infection Status",
       y = "Read Count")
print(plot2)

plot3 <- ggplot(infected_counts, aes(x = cell_number, y = z_score)) +
  geom_point() +
  geom_hline(yintercept = 3, linetype = "
  
plot1 <- ggplot(read_counts, aes(x = count, fill = infected)) +
  geom_histogram(binwidth = 50, position = "dodge", alpha = 0.7) +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Read Counts Distribution",
       x = "Read Count",
       y = "Frequency",
       fill = "Infection Status")
print(plot1)

plot2 <- ggplot(read_counts, aes(x = infected, y = count, fill = infected)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Boxplot of Read Counts by Infection Status",
       x = "Infection Status",
       y = "Read Count")
print(plot2)

plot3 <- ggplot(infected_counts, aes(x = cell_number, y = z_score)) +
  geom_point() +
  geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Z-Scores of Infected Read Counts",
       x = "Cell Number",
       y = "Z-Score")
print(plot3)


