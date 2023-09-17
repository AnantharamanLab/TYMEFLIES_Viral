# Read the data from the provided document
data <- read.table("UBA10906_virus_and_MAG_plots/pmoC_date2valid_pair.for_UBA10906.txt", header = FALSE, sep = "\t")

# Calculate the ratios AMG-containing virus/host and no-AMG-containing virus/host 
data$AMG_virus_host_Ratio <- data$V3 / data$V2
data$no_AMG_virus_host_Ratio <- data$V4 / data$V2

# Combine both ratios into a single data frame
combined_data <- data.frame(Ratio_Type = rep(c("AMG-containing virus/host", "no-AMG-containing virus/host"), each = nrow(data)),
                            Ratio_Value = c(data$AMG_virus_host_Ratio, data$no_AMG_virus_host_Ratio))

# Load the required libraries for plotting
library(ggplot2)
library(scales)

# Calculate mean values for each ratio type
mean_values <- aggregate(Ratio_Value ~ Ratio_Type, combined_data, mean)

# Generate a single boxplot for both AMG-containing virus/host and no-AMG-containing virus/host ratios
plot <- ggplot(combined_data, aes(x = Ratio_Type, y = Ratio_Value, fill = Ratio_Type)) +
  geom_boxplot(alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 1.5, shape = 16) +
  geom_text(data = mean_values, aes(x = Ratio_Type, y = Ratio_Value, label = sprintf("%.2f", Ratio_Value)), 
            vjust = -0.5, size = 3) +
  scale_y_continuous(trans = "log10", labels = comma) +
  labs(x = NULL, y = "Ratio Value (comma-separated)") +
  ggtitle("Boxplots for AMG-containing virus/host and no-AMG-containing virus/host Ratios") +
  theme_minimal() +
  scale_fill_manual(values = c("#0099CC", "#336666"))

# Save the plot as a PDF with height 5 and width 3 inches
ggsave("UBA10906_virus_and_MAG_plots/virus2host_ratio_boxplots.pdf", plot, height = 5, width = 3.5, units = "in")
