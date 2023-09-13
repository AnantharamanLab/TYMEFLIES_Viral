library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(plyr)

# Step 1 Plot the KO2tax2abun_fraction of 8 low diversity KOs
table <-read.table("AMG_analysis/KO2tax2abun_fraction.txt",sep = "\t",head=T)
table.subset <- table  %>% select(Head,K00036,K02294,K15895,K20811,K00033,K02706,K13017,K02703)
row.names(table.subset) <- table.subset$Head
table.subset[1] <- NULL
table.subset <- as.matrix(table.subset)
table.subset <-prop.table(table.subset,2) # Convert to percentage according to column sum
table.subset <- as.data.frame(table.subset)
table.subset$Head <- rownames(table.subset)

table.subset$Head <- revalue(table.subset$Head, c("NA;NA"="UNA;NA")) # Replace value

## Reshape the data
table.subset <- table.subset %>% gather("KO.ID", "Percentage",  1:(ncol(table.subset)-1))
## Get palette
getPalette <- colorRampPalette(brewer.pal(52, "Set1"))

p <- ggplot(table.subset, aes(fill=Head, y=Percentage, x=KO.ID)) + 
      geom_bar(position="stack", stat="identity")+
      # Make it pretty:
      theme_bw()+
      theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            legend.position="bottom")+
      # Use a nice color scheme:
      scale_fill_manual(values = colorRampPalette(brewer.pal(52, "Accent"))(52)) + xlab("KO ID")

p
ggsave(p,file="./AMG_analysis/KO2tax2abun_fraction.8_low_diversity_KOs.pdf", width = 11, height = 8.5, units = "in")


# Step 2 Plot the KO2host_tax2abun_fraction of 8 low diversity KOs
table2 <-read.table("AMG_analysis/KO2host_tax2abun_fraction.txt",sep = "\t",head=T)
table2.subset <- table2  %>% select(Head,K00036,K02294,K15895,K20811,K00033,K02706,K13017,K02703)
row.names(table2.subset) <- table2.subset$Head
table2.subset[1] <- NULL
table2.subset <- as.matrix(table2.subset)
table2.subset <-prop.table(table2.subset,2) # Convert to percentage according to column sum
table2.subset <- as.data.frame(table2.subset)
table2.subset$Head <- rownames(table2.subset)
table2.subset$Head <- revalue(table2.subset$Head, c("o__;f__"="Uo__;f__")) # Replace value

write.table(table2.subset,"./AMG_analysis/KO2host_tax2abun_fraction.8_low_diversity_KOs.txt",sep = "\t", quote=F);


# Read the modified table2
table2.mdf <-read.table("AMG_analysis/KO2host_tax2abun_fraction.8_low_diversity_KOs.mdf.txt",sep = "\t",head=T)
summary(table2.mdf) # right now your KO ID are columns, you need to gather them first:
table2.mdf <- table2.mdf %>% gather(KO.ID, Percentage, 2:ncol(table2.mdf))

# Define the desired order of KO IDs
desired_order <- c("K00033", "K00036", "K02703", "K02706", "K02294", "K20811", "K15895", "K13017")

# Create a factor variable with the desired order
table2.mdf$KO.ID <- factor(table2.mdf$KO.ID, levels = desired_order)

# Rest of your code
getPalette <- brewer.pal(10, "Set3")

p2 <- ggplot(table2.mdf, aes(fill = Head, y = Percentage, x = KO.ID)) + 
  geom_bar(position = "stack", stat = "identity") +
  # Make it pretty:
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom") +
  # Use a nice color scheme:
  scale_fill_manual(values = getPalette) + xlab("KO ID")

ggsave(p2, file = "./AMG_analysis/KO2host_tax2abun_fraction.8_low_diversity_KOs.pdf", width = 11, height = 8.5, units = "in")

