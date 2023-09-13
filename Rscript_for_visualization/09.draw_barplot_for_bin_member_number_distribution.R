library(tidyverse)
library(reshape2)
library(ggpubr)

# Read table
table <- read.table ("Bin_member_number_distribution.txt", sep = "\t", head = F, na.strings = c("NA"))
head(table)

# Add name to columns
colnames(table)[1] <- "Bin_member_number"
colnames(table)[2] <- "Frequency"

# Plot boxplot
p<-ggplot(table, aes(x=Bin_member_number, y=Frequency)) +
   geom_bar(stat="identity")+
   labs(title="Bin member number frequency",x="Bin member number", y = "Frequency (%)")+
   scale_x_discrete(limits=c("Bin2","Bin3","Bin4","Bin5","Bin6","Bin7","Bin8","Bin9","Bin10","Bin11","Bin12","Bin13","Bin14"))
   # Only keep Bin14 and before values that are > 1%
p

ggsave(p,file=paste0("Bin_member_number_distribution.pdf"), width = 12, height = 15, units = "cm")
ggsave(p,file=paste0("Bin_member_number_distribution.png"), width = 12, height = 15, units = "cm")
