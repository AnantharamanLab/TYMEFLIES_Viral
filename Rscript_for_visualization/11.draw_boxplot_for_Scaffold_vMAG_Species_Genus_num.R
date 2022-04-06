library(tidyverse)
library(reshape2)
library(ggpubr)

# Read table
table <- read.table ("Scaffold_vMAG_Species_Genus_num.txt", sep = "\t", head = F, na.strings = c("NA"))
head(table)

# Add name to columns
colnames(table)[1] <- "Item"
colnames(table)[2] <- "Number"

# Plot stack boxplot
p<-ggplot(table, aes(x=Item, y=Number)) +
  geom_bar(stat="identity")+
  labs(title="Item number",x="Item", y = "Number")+
  scale_x_discrete(limits=c("Viral scaffold number","vMAG number","Species-level vOTU number","Genus-level vOTU number"))+
  scale_fill_discrete(limits = c("Viral scaffold number","vMAG number","Species-level vOTU number","Genus-level vOTU number"))
p

ggsave(p,file=paste0("Scaffold_vMAG_Species_Genus_num.pdf"), width = 10, height = 10, units = "cm")
ggsave(p,file=paste0("Scaffold_vMAG_Species_Genus_num.png"), width = 10, height = 10, units = "cm")
