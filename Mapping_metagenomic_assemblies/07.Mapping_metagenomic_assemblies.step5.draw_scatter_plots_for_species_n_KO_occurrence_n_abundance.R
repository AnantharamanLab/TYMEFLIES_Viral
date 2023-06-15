library(ggplot2)
library(ggpmisc)

# Step 1 For scatter plot of Species2occurrence_n_abundance
## Read table
table <- read.table ("./AMG_analysis/Species2occurrence_n_abundance.txt", sep = "\t", head = F, na.strings = c("NA"))
head(table)

## Make scatter plot
p <- ggplot(table, aes(x=V2, y=V3)) + 
     geom_point(color="#f8766d", alpha = 0.6, size = 0.35)+
     stat_smooth(method = lm, color = "grey", size = 0.35)
p
ggsave(p,file="./AMG_analysis/Species2occurrence_n_abundance.pdf", width = 9, height = 7, units = "cm")

# Step 2 For scatter plot of KO2occurrence_n_abundance
## Read table
table2 <- read.table ("./AMG_analysis/KO2occurrence_n_abundance.txt", sep = "\t", head = F, na.strings = c("NA"))
head(table2)

## Make scatter plot
p2 <- ggplot(table2, aes(x=V2, y=V3)) + 
      geom_point(color="#f8766d", alpha = 0.6, size = 0.35)+
      stat_smooth(method = lm, color = "grey", size = 0.35)+
      stat_poly_eq(parse=T, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), formula=y~x)+
      geom_text(label=table2$V1, hjust=0.5, vjust=0.3, size = 0.3)
p2
ggsave(p2,file="./AMG_analysis/KO2occurrence_n_abundance.pdf", width = 9, height = 7, units = "cm")
