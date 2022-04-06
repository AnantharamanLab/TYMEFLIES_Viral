library(tidyverse)
library(reshape2)
library(ggpubr)

# Read table
table <- read.table ("Sample_num2species_num_for_refraction_curve.txt", sep = "\t", head = T, na.strings = c("NA"))
head(table)

# Make line plot
p<- ggplot(table, aes(x=Sample, y=Mean, color = "brown")) + 
    geom_line() +
    geom_point(size = 0.5, alpha = 0.25)+
    geom_errorbar(aes(ymin=Mean-Std, ymax=Mean+Std), width=.2, alpha = 0.25,
                position=position_dodge(0.05))
p

ggsave(p,file=paste0("Sample_num2species_num_for_refraction_curve.pdf"), width = 15, height = 12, units = "in")
ggsave(p,file=paste0("Sample_num2species_num_for_refraction_curve.png"), width = 15, height = 12, units = "in")