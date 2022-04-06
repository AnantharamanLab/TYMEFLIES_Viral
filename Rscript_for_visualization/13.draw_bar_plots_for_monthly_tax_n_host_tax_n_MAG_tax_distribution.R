library(tidyverse)
library(reshape2)
library(ggpubr)
library(RColorBrewer)

# Step 1 Plot Month2family2abun bar plots
# Read table
table <- read.table ("Month2family2abun.mdfed.txt", sep = "\t", head = T, na.strings = c("NA"),row.names = 1)
head(table)
table <- as.matrix(table)
table <-prop.table(table,2) # Convert to percentage according to column sum
table <- as.data.frame(table)
table$Head <- rownames(table)

# Write percentage table
write.table(table,"Month2family2abun.percentage.txt",sep = "\t",quote=F)

# Read percentage table (mdfed)
# Remove families that always < 1% across all months and add one column "XMean"
table <- read.table ("Month2family2abun.percentage.mdfed.txt", sep = "\t", head = T, na.strings = c("NA"))

## Reshape the data
table <- table %>% gather("Month", "Percentage",  2:ncol(table))
## Get palette
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

p <- ggplot(table, aes(fill=Head, y=Percentage, x=Month)) + 
  geom_bar(position="stack", stat="identity")+
  # Make it pretty:
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom",
        text = element_text(size = 5))+
  # Use a nice color scheme:
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Accent"))(9)) + xlab("Month")

p
ggsave(p,file="./Month2family2abun.pdf", width = 11, height = 8.5, units = "cm")

# Step 2 Plot Month2host_family2abun bar plots
# Read table
table2 <- read.table ("Month2host_family2abun.mdfed.txt", sep = "\t", head = T, na.strings = c("NA"),row.names = 1)
head(table2)
table2 <- as.matrix(table2)
table2 <-prop.table(table2,2) # Convert to percentage according to column sum
table2 <- as.data.frame(table2)
table2$Head <- rownames(table2)

# Write percentage table
write.table(table2,"Month2host_family2abun.percentage.txt",sep = "\t",quote=F)

# Read percentage table (mdfed)
# Remove families that always < 0.3% across all months and add one column "XMean"
table2 <- read.table ("Month2host_family2abun.percentage.mdfed.txt", sep = "\t", head = T, na.strings = c("NA"))

## Reshape the data
table2 <- table2 %>% gather("Month", "Percentage",  2:ncol(table2))
## Get palette
getPalette <- colorRampPalette(brewer.pal(23, "Set1"))

p2 <- ggplot(table2, aes(fill=Head, y=Percentage, x=Month)) + 
  geom_bar(position="stack", stat="identity")+
  # Make it pretty:
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom",
        text = element_text(size = 3))+
  # Use a nice color scheme:
  scale_fill_manual(values = colorRampPalette(brewer.pal(23, "Accent"))(23)) + xlab("Month")

p2
ggsave(p2,file="./Month2host_family2abun.pdf", width = 11, height = 8.5, units = "cm")

# Step 3 Plot Month2family2abun (for MAGs) bar plots
# Read table
table3 <- read.table ("../MAG_abundance/Family2month2abun_for_MAG.mdfed.txt", sep = "\t", head = T, na.strings = c("NA"),row.names = 1)
head(table3)
table3 <- as.matrix(table3)
table3 <-prop.table(table3,2) # Convert to percentage according to column sum
table3 <- as.data.frame(table3)
table3$Head <- rownames(table3)

# Write percentage table
write.table(table3,"Month2family2abun_for_MAG.percentage.txt",sep = "\t",quote=F)

# Read percentage table (mdfed)
# Remove families that always < 2.5% across all months and add one column "XMean"
table3 <- read.table ("Month2family2abun_for_MAG.percentage.mdfed.txt", sep = "\t", head = T, na.strings = c("NA"))

## Reshape the data
table3 <- table3 %>% gather("Month", "Percentage",  2:ncol(table3))
## Get palette
getPalette <- colorRampPalette(brewer.pal(41, "Set1"))

p3 <- ggplot(table3, aes(fill=Head, y=Percentage, x=Month)) + 
  geom_bar(position="stack", stat="identity")+
  # Make it pretty:
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom",
        text = element_text(size = 3))+
  # Use a nice color scheme:
  scale_fill_manual(values = colorRampPalette(brewer.pal(41, "Accent"))(41)) + xlab("Month")

p3
ggsave(p3,file="./Month2family2abun_for_MAG.pdf", width = 11, height = 11, units = "cm")
