library(tidyverse)
library(reshape2)
library(ggpubr)

# Read table
table <- read.table ("Bin2checkv_quality_and_length.txt", sep = "\t", head = F, na.strings = c("NA"))
head(table)

# Add name to columns
colnames(table)[1] <- "Phage genomes"
colnames(table)[2] <- "Checkv_quality"
colnames(table)[3] <- "Length"

# Plot boxplot

p <- ggplot(data=table,aes(x=Checkv_quality,y=Length))+
  geom_boxplot(aes(fill=Checkv_quality),outlier.shape = NA)+
  scale_y_log10()+
  labs(title="CheckV quality to phage genome length",x="CheckV quality", y = "Length (bp)")+
  #geom_jitter(width = 0.25, alpha = 0.1, color = 'black')+ # Too many dots to plot (cause too big picture)
  scale_x_discrete(limits=c("Complete","High-quality","Medium-quality","Low-quality","Not-determined"))+
  scale_fill_discrete(name = "Group", limits = c("Complete","High-quality","Medium-quality","Low-quality","Not-determined"))+
  stat_summary(fun=mean, colour="darkred", geom="point", hape=18, size=3,show.legend = FALSE)+
  stat_summary(fun=mean, colour="black", geom="text", show.legend = FALSE, 
               vjust=-1, aes( label=round(10^..y.., digits=1)))

table.summary <- tapply(table$Length, table$Checkv_quality, summary)
table.summary

ggsave(p,file=paste0("Bin2checkv_quality_and_length.pdf"), width = 20, height = 15, units = "cm")
ggsave(p,file=paste0("Bin2checkv_quality_and_length.png"), width = 20, height = 15, units = "cm")