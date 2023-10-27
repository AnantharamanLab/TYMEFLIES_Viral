# Load required libraries
library(ggplot2)

# Given data
genes <- c("psbA", "psbD", "pmoC", "katG", "ahbD", "gpx", "cobS", "cobT", "nadE", "cysC", "cysD", "cysH", "cysK")
percentages <- data.frame(
  Genes = factor(genes, levels = genes),
  `A75-100` = c(91.5, 94.2, 89.3, 86.7, 85.4, 89.3, 85.9, 89.8, 100, 91.7, 98.3, 92.7, 92.7),
  `B50-75` = c(82.6, 81.6, 58.3, 63.6, 82.3, 80.2, 77.8, 78.7, 91.5, 89.4, 68.4, 84.7, 82.9)
)

# Convert data to long format for ggplot2
percentages_long <- tidyr::gather(percentages, Range, Percentage, -Genes)

# Create the bar plot
bar_plot <- ggplot(percentages_long, aes(x = Genes, y = Percentage, fill = Range)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.85, color = "white") +
  scale_fill_manual(values = c("#20B2AA", "#ADD8E6")) +
  labs(title = "Gene Expression Percentages",
       x = "Genes",
       y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill = "white", size = 0.5),
        legend.box.background = element_rect(color = "black", size = 0.5))

# Print the plot
print(bar_plot)
