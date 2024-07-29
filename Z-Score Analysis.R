#codes incorporating the 0 value 

library(MASS)
library(ggplot2)
library(pscl)  # For zero-inflated models as there were some 0 value

# Sample data
read_counts <- data.frame(
  cell_number = 1:96,
  count = c(
0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 2, 1, 1, 3, 1, 1, 0, 5, 1, 0, 0, 1, 0, 2, 2, 0, 0, 4, 0, 1, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 1, 0, 10, 1, 0, 5872, 0, 1, 0, 2, 0, 1, 1, 0, 0, 75, 16, 4099, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 5, 0, 3, 0, 0, 0, 0
)
)

read_counts$infected <- ifelse(read_counts$cell_number %in% c(1:24, 73:96), "Uninfected", "Infected")

uninfected_counts <- read_counts[read_counts$infected == "Uninfected", ]
infected_counts <- read_counts[read_counts$infected == "Infected", ]

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

# Mann-Whitney U Test with the threshold of 1000 reads 
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
