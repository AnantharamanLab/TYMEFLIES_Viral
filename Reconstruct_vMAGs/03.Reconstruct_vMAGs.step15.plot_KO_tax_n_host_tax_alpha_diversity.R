library(ggplot2)
library(tidyverse)
library(vegan)
library(ggpmisc)

# Step 1 Plot the scatter plot of KO alpha diversity (based on family) to KO occurrence
KO.family <-read.table("AMG_analysis/KO2family2viral_gn_num.txt",sep = "\t",head = T, row.names = 1)
KO.family.shannon <- diversity(KO.family, index = "shannon")
KO.family.shannon <- as.data.frame(KO.family.shannon)
KO.family.shannon$KO <- rownames(KO.family.shannon)

KO.family.simpson <- diversity(KO.family, index = "simpson")
KO.family.simpson <- as.data.frame(KO.family.simpson)
KO.family.simpson$KO <- rownames(KO.family.simpson)

## Store occurrence 
KO.occurrence_n_abundance <- read.table ("AMG_analysis/KO2occurrence_n_abundance.txt",sep = "\t",head = F)
colnames(KO.occurrence_n_abundance) <- c('KO','Occurrence','Abundance')

## Merge two tables
KO.family.shannon2occurrence <- left_join(KO.family.shannon, KO.occurrence_n_abundance, by = "KO")   

p <- ggplot(KO.family.shannon2occurrence, aes(x=Occurrence, y=KO.family.shannon)) + 
     geom_point(color="#f8766d", alpha = 0.6, size = 0.35)
p
ggsave(p,file="AMG_analysis/KO.family.shannon2occurrence.pdf", width = 9, height = 7, units = "cm")

KO.family.simpson2occurrence <- left_join(KO.family.simpson, KO.occurrence_n_abundance, by = "KO")   

p2 <- ggplot(KO.family.simpson2occurrence, aes(x=Occurrence, y=KO.family.simpson)) + 
  geom_point(color="#f8766d", alpha = 0.6, size = 0.35)
p2
ggsave(p2,file="AMG_analysis/KO.family.simpson2occurrence.pdf", width = 9, height = 7, units = "cm")

# Step 2 Plot the scatter plot of KO alpha diversity (based on host family) to KO occurrence
KO.host_family <-read.table("AMG_analysis/KO2host_family2viral_gn_num.txt",sep = "\t",head = T, row.names = 1)

KO.host_family.shannon <- diversity(KO.host_family, index = "shannon")
KO.host_family.shannon <- as.data.frame(KO.host_family.shannon)
KO.host_family.shannon$KO <- rownames(KO.host_family.shannon)

KO.host_family.simpson <- diversity(KO.host_family, index = "simpson")
KO.host_family.simpson <- as.data.frame(KO.host_family.simpson)
KO.host_family.simpson$KO <- rownames(KO.host_family.simpson)

## Merge two tables
KO.host_family.shannon2occurrence <- left_join(KO.host_family.shannon, KO.occurrence_n_abundance, by = "KO")   

p3 <- ggplot(KO.host_family.shannon2occurrence, aes(x=Occurrence, y=KO.host_family.shannon)) + 
      geom_point(color="#f8766d", alpha = 0.6, size = 0.35)+
      stat_smooth(method = lm, color = "grey", size = 0.35)+
      geom_text(label=KO.host_family.shannon2occurrence$KO, hjust=0.5, vjust=0.3, size = 0.3)
p3
ggsave(p3,file="AMG_analysis/KO.host_family.shannon2occurrence.pdf", width = 9, height = 7, units = "cm")

KO.host_family.simpson2occurrence <- left_join(KO.host_family.simpson, KO.occurrence_n_abundance, by = "KO")   

p4 <- ggplot(KO.host_family.simpson2occurrence, aes(x=Occurrence, y=KO.host_family.simpson)) + 
  geom_point(color="#f8766d", alpha = 0.6, size = 0.35)+
  stat_smooth(method = lm, color = "grey", size = 0.35)+
  stat_poly_eq(parse=T, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), formula=y~x)+
  geom_text(label=KO.host_family.simpson2occurrence$KO, hjust=0.5, vjust=0.3, size = 0.3)

p4
ggsave(p4,file="AMG_analysis/KO.host_family.simpson2occurrence.pdf", width = 9, height = 7, units = "cm")
