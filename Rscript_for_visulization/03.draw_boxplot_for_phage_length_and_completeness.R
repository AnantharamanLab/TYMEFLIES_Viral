library(tidyverse)
library(reshape2)
library(ggpubr)
install.packages("purrr")
library("purrr")

# Read table
table <- read.table ("phage_scaffold_and_vMAGs_length_and_completeness_info.txt", sep = "\t", head = F, na.strings = c("NA"))
head(table)

# Add name to columns
colnames(table)[1] <- "Phage scaffolds or bins"
colnames(table)[2] <- "Group"
colnames(table)[3] <- "Length"
colnames(table)[4] <- "Completeness"

# Plot boxplot
## For length
table.length <- table[,c(1,2,3)]

## Pick compared groups
compaired <- list(c("Phage scaffolds", "vMAGs + unbinned scaffolds"),c("Phage scaffolds", "vMAGs"))

p.length <- ggplot(data=table.length,aes(x=Group,y=Length))+
  geom_boxplot(aes(fill=Group),outlier.shape = NA)+
  scale_y_log10()+
  labs(title="Phage scaffold or bin length",x="Group", y = "Length (bp)")+
  #geom_jitter(width = 0.25, alpha = 0.1, color = 'black')+ # Too many dots to plot (cause too big picture)
  scale_x_discrete(limits=c("Phage scaffolds","vMAGs + unbinned scaffolds","vMAGs","Binned scaffolds (within vMAGs)","Unbinned scaffolds"))+
  scale_fill_discrete(name = "Group", limits = c("Phage scaffolds","vMAGs + unbinned scaffolds","vMAGs","Binned scaffolds (within vMAGs)","Unbinned scaffolds"))+
  stat_summary(fun=mean, colour="darkred", geom="point", hape=18, size=3,show.legend = FALSE)+
  stat_summary(fun=mean, colour="red", geom="text", show.legend = FALSE, 
                 vjust=-0.7, aes( label=round(10^..y.., digits=1)))+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)

table.length.summary <- tapply(table.length$Length, table.length$Group, summary)
table.length.summary

## For completeness
table.completeness <- table[,c(1,2,4)]

p.completeness <- ggplot(data=table.completeness,aes(x=Group,y=Completeness))+
  geom_boxplot(aes(fill=Group),outlier.shape = NA)+
  labs(title="Phage scaffold or bin completeness",x="Group", y = "Completeness (%)")+
  #geom_jitter(width = 0.25, alpha = 0.1, color = 'black')+
  scale_x_discrete(limits=c("Phage scaffolds","vMAGs + unbinned scaffolds","vMAGs","Binned scaffolds (within vMAGs)","Unbinned scaffolds"))+
  scale_fill_discrete(name = "Group", limits = c("Phage scaffolds","vMAGs + unbinned scaffolds","vMAGs","Binned scaffolds (within vMAGs)","Unbinned scaffolds"))+
  stat_summary(fun=mean, colour="darkred", geom="point", hape=18, size=3,show.legend = FALSE)+
  stat_summary(fun=mean, colour="red", geom="text", show.legend = FALSE, 
               vjust=-0.7, aes( label=round(..y.., digits=1)))+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)

table.completeness.summary <- tapply(table.completeness$Completeness, table.length$Group, summary)
table.completeness.summary

ggarrange(p.length,p.completeness,labels=c("A","B"), ncol=2, nrow=1, widths = c(1,1))

ggsave(file=paste0("Phage_scaffold_or_bin_length_and_completeness.pdf"), width = 55, height = 15, units = "cm")
ggsave(file=paste0("Phage_scaffold_or_bin_length_and_completeness.png"), width = 55, height = 15, units = "cm")

