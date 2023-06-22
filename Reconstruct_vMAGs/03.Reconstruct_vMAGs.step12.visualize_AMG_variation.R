library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

# Step 1 Plot the general statistics of vOTU and ko combinations
# Note that this "Species_level_vOTU_AMG_variation_statistics.table_1.txt" is 
# originated from "Species_level_vOTU_AMG_variation_statistics.txt"
table <-read.table("AMG_analysis/Species_level_vOTU_AMG_variation_statistics.table_1.txt",sep = "\t",head=T)

p <- ggplot(table, aes(x = vOTU.size, y = Fraction, fill = AMG.KO.frequency.category, label = Fraction)) +
     geom_bar(stat = "identity") +
     geom_text(size = 3, position = position_stack(vjust = 0.5))+
     # Make it pretty:
     theme_bw()+
     theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position="bottom")
  
ggsave(p,file="./AMG_analysis/Species_level_vOTU_AMG_variation_statistics.table_1.pdf", width = 10, height = 16, units = "in")

# Step 2 Plot the scatter plot for KO2ko_abun_n_mean_ko_freq.txt
## Read table (In the mdfed table, HighOccurrence information was added - occurrence higher than 400 )
table2 <- read.table ("./AMG_analysis/KO2ko_abun_n_mean_ko_freq.mdfed.txt", sep = "\t", head = F, na.strings = c("NA"))
head(table2)

## Make scatter plot
p2 <- ggplot(table2, aes(x=V3, y=V2, color = V4)) + 
      geom_point(alpha = 0.6, size = 0.35)+
      geom_text(label=table2$V1, hjust=0.5, vjust=0.3, size = 0.3)
p2
p2 <- p2 + theme(legend.position = "none")
ggsave(p2,file="./AMG_analysis/KO2ko_abun_n_mean_ko_freq.pdf", width = 9, height = 7, units = "cm")

