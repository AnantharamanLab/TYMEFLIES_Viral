library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

# Step 1 Plot the general statistics of vOTU and ko combinations
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

# Step 2 Plot the KO2vOTU_tax2abun.txt result
table2 <-read.table("AMG_analysis/Species_level_vOTU_AMG_variation_statistics.table_2.txt",sep = "\t",head=T)
## Reshape the data
table2 <- table2 %>% gather("KO.name", "Percentage", 2:ncol(table2))
## Get palette
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

p2 <- ggplot(table2, aes(fill=Head, y=Percentage, x=KO.name)) + 
      geom_bar(position="stack", stat="identity")+
      # Make it pretty:
      theme_bw()+
      theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            legend.position="bottom")+
      # Use a nice color scheme:
      scale_fill_manual(values = colorRampPalette(brewer.pal(18, "Accent"))(18)) + xlab("KO Name")

ggsave(p2,file="./AMG_analysis/Species_level_vOTU_AMG_variation_statistics.table_2.pdf", width = 11, height = 8.5, units = "in")


# Step 3 Plot the KO2vOTU_host_tax2abun.txt result
table3 <-read.table("AMG_analysis/Species_level_vOTU_AMG_variation_statistics.table_3.txt",sep = "\t",head=T)
## Reshape the data
table3 <- table3 %>% gather("KO.name", "Percentage", 2:ncol(table3))
## Get palette
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

p3 <- ggplot(table3, aes(fill=Head, y=Percentage, x=KO.name)) + 
  geom_bar(position="stack", stat="identity")+
  # Make it pretty:
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom")+
  # Use a nice color scheme:
  scale_fill_manual(values = colorRampPalette(brewer.pal(18, "Accent"))(18)) + xlab("KO Name")

ggsave(p3,file="./AMG_analysis/Species_level_vOTU_AMG_variation_statistics.table_3.pdf", width = 11, height = 8.5, units = "in")

# Step 4 Plot the scatter plot for KO2ko_abun_n_mean_ko_freq.txt
## Read table
table4 <- read.table ("./AMG_analysis/KO2ko_abun_n_mean_ko_freq.mdfed.txt", sep = "\t", head = F, na.strings = c("NA"))
head(table4)

## Make scatter plot
p4 <- ggplot(table4, aes(x=V3, y=V2, color = V4)) + 
      geom_point(alpha = 0.6, size = 0.35)+
      geom_text(label=table4$V1, hjust=0.5, vjust=0.3, size = 0.3)
p4
p4 <- p4 + theme(legend.position = "none")
ggsave(p4,file="./AMG_analysis/KO2ko_abun_n_mean_ko_freq.pdf", width = 9, height = 7, units = "cm")

# Step 5 Plot the seasonal KO-carrying viral genome distribution by "KO2dates_in_a_year.mdf.txt"
## Read table
table5 <- read.table ("./AMG_analysis/KO2dates_in_a_year.mdfed.txt", sep = "\t", head = T, na.strings = c("NA"))
head(table5)

## Only keep high occurrence KOs
table5.subset <- subset(table5, OccurrenceState %in% c("HighOccurrenceKO"))

## Only keep AMG KO fraction >= 1.25%
table5.subset <- left_join(table5.subset, table4, by = c("Head" = "V1"))
table5.subset["1", "V2"] <- 1
table5.subset["2", "V2"] <- 1
table5.subset <- table5.subset %>% filter(V2 >= 0.0125)
row.names(table5.subset) <- table5.subset$Head
table5.subset = select(table5.subset, -1:-2, -368:-370)
table5.subset <- as.matrix(table5.subset)

## Make heatmap
p5 <- pheatmap(table5.subset,cluster_row = F, cluster_col= F,fontsize = 3,color = colorRampPalette(c("White","#f8766d"))(2),cellwidth = 0.4, cellheight = 5)

ggsave(p5,file="./AMG_analysis/KO2dates_in_a_year.pdf", width = 11, height = 8.5, units = "cm")