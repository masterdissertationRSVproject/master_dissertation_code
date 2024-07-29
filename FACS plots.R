#Codes to visualise FACS result based on the raw data
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("flowCore")

install.packages("ggplot2")

library(flowCore)
library(ggplot2)

fcs_directory <- "/path_to_my_files/RSV Gene Expression/FACS"

fcs_files <- list.files(path = fcs_directory, pattern = "\\.fcs$", full.names = TRUE)

all_fcs_data <- list()

# Loop through each FCS file
for (fcs_file in fcs_files) {

  fcs_data <- read.FCS(fcs_file)
  
  data <- exprs(fcs_data)
  
  data <- as.data.frame(data)
  data$file_id <- basename(fcs_file)
  
  all_fcs_data[[basename(fcs_file)]] <- data
  
  cat("Processed:", basename(fcs_file), "\n")
}

# Combine all data into a single data frame
combined_fcs_data <- do.call(rbind, all_fcs_data)

head(combined_fcs_data)

# Dot plot FSC vs SSC for all files 
ggplot(combined_fcs_data, aes(x = `FSC-A`, y = `SSC-A`, color = file_id)) +
  geom_point(alpha = 0.5) +
  labs(title = "FSC vs SSC Plot for All Files",
       x = "Forward Scatter (FSC-A)",
       y = "Side Scatter (SSC-A)") +
  theme_minimal()

for (file_id in unique(combined_fcs_data$file_id)) {
  file_data <- subset(combined_fcs_data, file_id == file_id)
  
  plot <- ggplot(file_data, aes(x = `FSC-A`, y = `SSC-A`)) +
    geom_point(alpha = 0.5) +
    labs(title = paste("FSC vs SSC Plot for", file_id),
         x = "Forward Scatter (FSC-A)",
         y = "Side Scatter (SSC-A)") +
    theme_minimal()
  
  ggsave(filename = paste0("FSC_vs_SSC_", file_id, ".png"), plot = plot)
}
