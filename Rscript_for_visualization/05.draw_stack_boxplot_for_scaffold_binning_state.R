library(tidyverse)
library(reshape2)
library(ggpubr)

# Read table
table <- read.table ("Scaffold_to_binning_state_before_and_after_binning.txt", sep = "\t", head = F, na.strings = c("NA"))
head(table)

# Add name to columns
colnames(table)[1] <- "Before_after_binning"
colnames(table)[2] <- "Binned_unbinned"
colnames(table)[3] <- "Scaffold_percentage"

# Plot stack boxplot
p <- ggplot(data=table, aes(x=Before_after_binning, y=Scaffold_percentage, fill=Binned_unbinned)) +
     geom_bar(stat="identity")+
     scale_x_discrete(limits=c("Before binning","After binning"))+
     scale_fill_discrete(name = "Binned/unbinned", limits = c("unbinned","binned"))+
     labs(title="Scaffold binned/unbinned percentage", x = "",y = "Scaffold percentage (%)")
p

ggsave(p,file=paste0("Scaffold_binned_unbinned_percentage.pdf"), width = 10, height = 10, units = "cm")
ggsave(p,file=paste0("Scaffold_binned_unbinned_percentage.png"), width = 10, height = 10, units = "cm")
