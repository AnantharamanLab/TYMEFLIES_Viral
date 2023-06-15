library(ggplot2)
library(reshape2)
library(scales)

library(ggplot2)
library(reshape2)
library(scales)

# Read the table data
table_data <- read.table("AMG_analysis/KO2season2KO_carrying_genome_distribution_percentage.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract KO numbers
ko_numbers <- table_data$KO

# Remove the first column (KO numbers) and convert to matrix
table_matrix <- as.matrix(table_data[, -1])

# Calculate the column averages
column_averages <- colMeans(table_matrix)

# Add the column averages as a new row to the matrix
table_matrix <- rbind(table_matrix, column_averages)

# Rename the row names of the matrix to match KO numbers
rownames(table_matrix) <- c(ko_numbers, "Average")

# Reshape the data into long format for ggplot
table_data_long <- melt(table_matrix)

# Set up the color scale
color_scale <- scale_fill_gradient(low = "orange", high = "red")

# Generate the heatmap using ggplot2
heatmap_plot <- ggplot(table_data_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  color_scale +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = "Heatmap") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(size = 12, face = "bold")) +
  geom_text(aes(label = percent(value, accuracy = 1)), color = "black", size = 3)

# Save the heatmap as a PDF file
ggsave("AMG_analysis/KO2season2KO_carrying_genome_distribution_percentage.pdf", heatmap_plot, width = 8, height = 6)
